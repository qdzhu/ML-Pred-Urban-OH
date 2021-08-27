function [ products, product_coeff, reactants, reactant_coeff, rate_fxn, isphoto ] = read_wrf_mech_line( tline_in )
%[ PRODUCTS, PRODUCT_COEFF, REACTANTS, REACTANT_COEFF, RATE_FXN ] = READ_WRF_MECH_LINE( TLINE_IN )
%   This function takes a string that is a line from the WRF-Chem KPP .eqn
%   file and separate it into products and their coefficients, reactants
%   and their coefficients, and the rate coefficient.
%
%   PRODUCTS and REACTANTS are cell arrays of the unique products and
%   reactants detected in the reaction. Their coefficients are returned as
%   vectors the same size as the cell array.  If a chemical species shows
%   up more than once on one side of the reaction, it will only occur once
%   in PRODUCTS or REACTANTS but its coefficient will reflect how many
%   times it appeared; e.g. H2O2 --> HO + HO will return HO as the only
%   product but with a coefficient of 2. Currently, if a third body is
%   specified as M{O2} or something similar, the part in {} will be ignored
%   and it will just be returned as an M.
%
%   RATE_FXN is a function handle to the appropriate rate constant
%   function. Non-photolytic rate constants will required the input of
%   temperature (in K) and number density of air (in molec./cm^3).
%
%   The KPP is described in Damian, et al. 2002; Computers and Chemical
%   Engineering, 26, 1567, http://www2.mmm.ucar.edu/wrf/WG2/GPU/KPP-sdarticle.pdf


E = JLLErrors;

% Get rid of the equation number which is enclosed in {}
i=strfind(tline_in,'}');
tline=tline_in(i+1:end);

% Split into equation and rate constant
tmp = strtrim(strsplit(tline,':'));
rate_const_in = tmp{2};
% Anything in curly braces is a comment and so can be removed. See Damian
% et al. section 4
tmp{1} = regexprep(tmp{1},'{.*}','');
% Then split into products and reactants
tmp = strtrim(strsplit(tmp{1},'='));
products = strtrim(strsplit(tmp{2},'+'));
reactants = strtrim(strsplit(tmp{1},'+'));
reactants = reactants(~strcmp(reactants,'hv'));

% Handle coefficients in reactants and products, or a product being
% specified multiple times.
[products, product_coeff] = process_eqn(products);
[reactants, reactant_coeff] = process_eqn(reactants);

% Parse the rate constant
% Get rid of anything from a semicolon to the end
sc = strfind(rate_const_in,';');
if ~isempty(sc)
    rate_const = rate_const_in(1:sc-1);
end
% Second remove any _dp's (which are used to indicate double precision),
% convert scientific notation of #.##D-## to #.##e-#, change calls to
% EXP to exp, and convert the fortran exponentiation operator (**) to the
% Matlab one (.^)
rate_const = regexprep(rate_const,'_dp','');
rate_const = regexprep(rate_const,'(?<=\d)D(?=[-\d])','e');
rate_const = regexprep(rate_const,'EXP','exp');
rate_const = regexprep(rate_const,'\*\*',' .^ ');

if strcmp(rate_const(1),'j')
    % Dealing with a photolysis rate constant. Just return a string like
    % 'Pj_no2' which will get converted into a function call later to
    % reference the output of TUV.
    [s,e] = regexp(rate_const, '\(.*\)');
    if ~isempty(s)
        rate_const = rate_const(s+1:e-1);
    end
    rate_fxn = rate_const;
    isphoto = true;
else
    % Get function handles to the possible rate functions
    ARR2 = wrf_rate_expr('ARR2');
    TROE = wrf_rate_expr('TROE');
    TROEE = wrf_rate_expr('TROEE');
    k46 = wrf_rate_expr('k46');
    ko1d = wrf_rate_expr('ko1d');
    rate_fxn = eval(sprintf('@(TEMP, C_M) %s',rate_const));
    isphoto = false;
end


end

function [chem_names_out, chem_coeff_out] = process_eqn(chem)
chem_names = chem;
chem_coeff = ones(size(chem_names));
for j=1:numel(chem_names)
    [s,e] = regexp(chem_names{j},'(?<!\w)\d*\.*\d*');
    if ~isempty(s)
        % Calculate the coefficient: if the species if mentioned more than
        % once, account for that.
        chem_coeff(j) = str2double(chem_names{j}(s:e));
        chem_names{j}(s:e)=[];
        chem_names{j} = strtrim(chem_names{j});
    end
end
chem_names_out = unique(chem_names);
chem_coeff_out = nan(size(chem_names_out));
for j=1:numel(chem_names_out)
    xx = strcmp(chem_names_out{j}, chem_names);
    chem_coeff_out(j) = sum(chem_coeff(xx));
end
end