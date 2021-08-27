function [ phot_rates ] = call_tuv( wrf_js, date_in, hour_in, lon_in, lat_in, utc_bool )
% [ PHOT_RATES ] = CALL_TUV( WRF_JS, DATE_IN, HOUR_IN, LON_IN, LAT_IN, UTC_BOOL )
%   Calls a compiled TUV model for the date and hour specified. DATE_IN
%   should be a valid date number or string, HOUR_IN should be an integer.
%   WRF_JS should be the list of photolysis rates needed.

E=JLLErrors;

if ~iscell(wrf_js)
    wrf_js = strsplit(wrf_js);
end

% Set up and call the TUV model
if ~ismac
    warning('System calls to TUV model may not work on non-mac computers as written')
end

if ismac
    % This was a kludge necessary to fix some error occurring with a
    % gfortran compiled TUV
    setenv('DYLD_LIBRARY_PATH', '/usr/local/bin:/opt/local/lib:')
end


mfile_dir = fileparts(mfilename('fullpath'));
tuv_dir = fullfile(mfile_dir, 'tuv5.2_source');


if utc_bool
    tmzone = 0;
else
    tmzone = round(lon_in/15);
end

set_date_in_tuv_input(date_in, hour_in, lon_in, lat_in, tmzone, tuv_dir);

wd = cd(tuv_dir);
status = system(fullfile(tuv_dir,'tuv'));
if status ~= 0
    E.callError('tuv_failed','TUV did not execute successfully');
end
cd(wd)

% Convert the wrf_j names into the TUV equations.
tuv_eqns = cell(size(wrf_js));
eqn_numbers = nan(size(wrf_js));
for a=1:numel(wrf_js)
    tuv_eqns{a} = wrf_tuv_eqn_mapper(wrf_js{a});
end

% Open the file output by TUV and read the resulting photolysis rates.
fid = fopen(fullfile(tuv_dir,'..','usrout.txt'));
tline = fgetl(fid);
read_rxns = false;
read_rates = false;
num_rxns = 0;
rate_table = [];
while ischar(tline)
    % are we done with the reactions section?
    if read_rxns && ~isempty(strfind(tline,'values at z'))
        read_rxns = false;
    end
    
    if read_rxns
        % These lines have the form '   # = rxn' so we read after the = and
        % compare to the WRF-TUV mapping, then record the reaction number.
        % Also keep track of how many reactions there are so we know how to
        % read the rate table.
        [s,e] = regexp(tline,'\d*(?= =)');
        num_rxns = str2double(tline(s:e));
        
        i = strfind(tline,'=');
        rxn = strtrim(tline(i+2:end));
        xx = ismember(tuv_eqns, rxn);
        if sum(xx) > 0
            eqn_numbers(xx) = num_rxns;
        end
    end
    
    if read_rates
        % import hour and SZA as floats, the rates in scientific notation
        format_spec = ['%f %f', repmat(' %e',1,num_rxns)];
        rate_line = sscanf(tline, format_spec);
        rate_table = cat(1, rate_table, rate_line');
    end
    
    % What section is next?
    if ~isempty(strfind(tline,'Photolysis rate coefficients'))
        read_rxns = true;
    elseif ~isempty(strfind(tline, 'time, hrs'))
        read_rates = true;
    end
    tline = fgetl(fid);
end
fclose(fid);

% Find the hour closest to the requested one and read the appropriate rates
[~,hh] = min(rate_table(:,1) - hour_in);
phot_rates = rate_table(hh, eqn_numbers+2); % the +2 accounts for the first two columns being hour and SZA

end

function set_date_in_tuv_input(date_in, hour_in, lon_in, lat_in, tmzone, tuv_dir)
fid = fopen(fullfile(tuv_dir,'INPUTS','usrinp-template'));
fidnew = fopen(fullfile(tuv_dir,'INPUTS','usrinp'),'w');
tline = fgetl(fid);
while ischar(tline)
    if ~isempty(strfind(tline,'tmzone'))
        newline = sprintf('lat = %14.3f   lon = %14.3f   tmzone = %11.1f',lat_in,lon_in,tmzone);
        fprintf(fidnew,'%s\n',newline);
    elseif ~isempty(strfind(tline,'iyear'))
        yr = sprintf('%12d',year(date_in));
        mn = sprintf('%11d',month(date_in));
        dy = sprintf('%13d',day(date_in));
        newline = sprintf('iyear = %s   imonth = %s   iday = %s',yr,mn,dy);
        fprintf(fidnew,'%s\n',newline);
    elseif ~isempty(strfind(tline,'tstart'))
        newline = sprintf('tstart = %1$11.3f   tstop = %1$12.3f   nt =               1\n',hour_in);
        fprintf(fidnew, newline);
    else
        fprintf(fidnew,'%s\n',tline);
    end
    tline = fgetl(fid);
end
fclose(fid);
fclose(fidnew);
end

