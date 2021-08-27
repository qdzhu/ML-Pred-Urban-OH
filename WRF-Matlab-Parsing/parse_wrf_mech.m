function [ J, species, isfixed, photo_calls ] = parse_wrf_mech( mech_name )
%PARSE_WRF_MECH Parses the specified WRF mechanism into Matlab
%   [ J, SPECIES, ISFIXED, PHOTO_CALLS ] = PARSE_WRF_MECH( 'r2smh' ) will
%   use the R2SMH mechanism developed by Azimeh Zare that builds off work
%   by Ellie Browne and other to provide a very comprehensive model of
%   alkyl nitrate NOx chemistry. Return variables are:
%
%       J: a sparse Jacobian essentially given as a cell array of anonymous
%       functions that each accept the inputs c, j, TEMP, and C_M which are
%       the concentration vector of species, the vector of photolysis rate
%       constants (see PHOTO_CALLS below), the temperature of the grid
%       cell, and the number density of air in the grid cell. These final
%       two should be scalar values.
%
%       SPECIES: a cell array of chemical names matching those in the
%       R2SMH/r2smh.spc file in the same order as given in J, that is, if
%       SPECIES{1} is 'NO' then J{1}(c, j, TEMP, C_M) will give the value
%       of d[NO]/dt.
%
%       ISFIXED: logical array that specifies which of the chemical species
%       has been explicitly set in the R2SMH mechanism to not change via
%       chemical reactions.  Usually these are species like H2O or general
%       third-body M which are controlled by meteorology, not chemistry. In
%       general, be sure to initialize a value for these.
%
%       PHOTO_CALLS: list of photolysis rate functions required, in the
%       order required. Pass to call_tuv to generate the vector j of
%       photolysis rates needed by the Jacobian.
%
%   [ ... ] = PARSE_WRF_MECH('r2smh-simple') uses a much simplified
%   mechanism focusing only on the main NOx cycle. Good for testing.
%
%   Josh Laughner <joshlaugh5@gmail.com> 27 June 2016

if ~exist('mech_name','var')
    mech_name = 'r2smh-simple';
end

mfile_dir = fileparts(mfilename('fullpath'));

species_file = fullfile(mfile_dir,'R2SMH',strcat(mech_name,'.spc'));
eqn_file = fullfile(mfile_dir,'R2SMH',strcat(mech_name,'.eqn'));

% The mechanism will be constructed as a vector of anonymous functions that
% accept a concentration vector and photolysis rate vector.
%
% Emissions must be handled separately in your model.

[species, isfixed] = parse_species(species_file);
[J, photo_calls] = construct_mechanism(eqn_file, species);
J(isfixed) = {@(c,j,TEMP,C_M) 0};
end

function [species, isfixed] = parse_species(species_file)
% This function will return the list of species defined for the mechanism
% as well as a logical vector indicating if those species are to be fixed
% by some physical process rather than the chemistry.
%
% First read in the species names from the species file, as well as if it
% is meant to be solved by chemistry or fixed by physical processes.
species = cell(1,1000); % will cut down in the end
isfixed = false(1,1000);
fixed_i = 0;
i=1;
fid = fopen(species_file);
tline = fgetl(fid);
while ischar(tline)
    tline = strtrim(tline);
    if ismember('#DEFVAR',tline) % the following species are variable
        fixed_i = 0;
    elseif ismember('#DEFFIX',tline) % the following species are fixed
        fixed_i = 1;
    elseif ismember('#',tline) % unanticipated definition bloc
        fprintf('Unknown definition block: %s\n',tline);
    else
        % Curly braces indicate a comment and should be removed. If the
        % line has only white space after that, skip it
        parsed_line = strtrim(regexprep(tline,'{.*}',''));
        if ~isempty(parsed_line)
            parsed_line = strrep(tline,'{','');
            parsed_line = strrep(parsed_line,'}','');
            parsed_line = strsplit(parsed_line,'='); % species comes before the = sign. ignore on the other side seems to just mean don't do mass balance checking (i.e. don't check that the number of atoms on both sides is conserved)
            spec = strtrim(parsed_line{1});
            xx = strcmp(spec, species);
            if sum(xx) > 0 && isfixed(xx) ~= fixed_i
                error('species_definition:species_fixed_and_var','Specie %s is defined as both a fixed and variable specie',spec)
            elseif sum(xx) > 0
                fprintf('Redundant definition of %s, skipping\n',spec);
            else
                species{i} = spec;
                isfixed(i) = fixed_i;
                i=i+1;
            end
        end
    end
    tline = fgetl(fid);
end
species(i:end) = [];
isfixed(i:end) = [];
fclose(fid);
end

function [J, photo_calls] = construct_mechanism(eqn_file, species)
J = cell(numel(species),1);
J(:) = {@(c,j,TEMP,C_M) 0};

fid = fopen(eqn_file);
tline = fgetl(fid);
photo_calls = {};
while ischar(tline)
    % skip comment lines or "dummy" reactions
    if ismember('#',tline) || ~isempty(regexp(tline,'JUNK','once'))
    else
        [products, product_coeff, reactants, reactant_coeff, k, isphoto] = read_wrf_mech_line(tline);
        
        % Construct the derivative function and add it to any existing such
        % functions in the mechanism for this species.  This relies heavily
        % on the fact that matlab function handles are fixed at the time of
        % creation, that is if you write:
        %   m = 2; b = 0;
        %   f = @(x) m.*x + b;
        % then f(1) will always return 2 even if you later change the value
        % of m or b.
        %
        % Crucially, this works for function handles too, so we can nest
        % old versions of handles in new ones, e.g.
        %   g = @() 1;
        %   g() % will return 1
        %   g = @() g() + 1;
        %   g() % will return 2
        %   g = @() g() + 1
        %   g() % will return 3
        rr = ismember(species, reactants);
        if sum(rr) ~= numel(reactants)
            nf = ~ismember(reactants, species);
            error('mechanism_parse:unknown_reactant','The reactant %s cannot be identified in the species list',strjoin(reactants(nf),', '));
        end
        
        % We need a function that accepts a consistent set of inputs
        % whether using a rate constant that is truly constant (i.e.
        % independent of temperature and number density of air), one that
        % depends on temperature and number density of air, or a photolysis
        % rate. f will be that function.
        if isphoto
            % We'll need this to know which function handles to construct
            % in the real mechanism.
            photo_calls{end+1} = k;
            k = sprintf('j(%d)',numel(photo_calls)); % This will call the proper element of vector j which must be passed to the jacobian and contains the photolysis rates. 
            f = @(c,j,TEMP,C_M) j(numel(photo_calls)) .* prod(c(rr) .^ reactant_coeff);
        else
            if ~isa(k,'function_handle')
                k = @(TEMP, C_M) k; % if k is a constant, make it a function that returns that constant but accepts inputs
            end
            f = @(c,j,TEMP,C_M) k(TEMP,C_M) .* prod(c(rr) .^ reactant_coeff);
        end
        

        
        % Add this for each reactant
        xx = find(ismember(species, reactants));
        if numel(xx) ~= numel(reactants)
            nf = ~ismember(reactants, species);
            error('mechanism_parse:unknown_reactant','The reactant %s cannot be identified in the species list',strjoin(reactants(nf),', '));
        end
        for i=1:numel(xx)
            RecursJ = J{xx(i)};
            rc_i = reactant_coeff(i);
            J{xx(i)} = @(c,j,TEMP,C_M) -rc_i .* f(c,j,TEMP,C_M) + RecursJ(c,j,TEMP,C_M);
        end
        
        % Add this for each product
        xx = find(ismember(species, products));
        if numel(xx) ~= numel(products)
            nf = ~ismember(products, species);
            error('mechanism_parse:unknown_reactant','The product %s cannot be identified in the species list',strjoin(products(nf),', '));
        end
        for i=1:numel(xx)
            RecursJ = J{xx(i)};
            pc_i = product_coeff(i);
            J{xx(i)} = @(c,j,TEMP,C_M) pc_i .* f(c,j,TEMP,C_M) + RecursJ(c,j,TEMP,C_M);
        end
    end
    tline = fgetl(fid);
end
fclose(fid);
end