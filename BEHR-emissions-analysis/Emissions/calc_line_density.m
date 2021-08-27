function [ no2_x, no2_linedens, no2_lindens_std, lon, lat, no2_mean, no2_std, num_valid_obs, nox, debug_cell ] = calc_line_density( fpath, fnames, center_lon, center_lat, theta, wind_logical, varargin )
%[ NO2_X, NO2_LINEDENS, NO2_LINEDENS_STD, LON, LAT, NO2_MEAN, NO2_STD, NUM_VALID_OBS] = CALC_LINE_DENSITY( FPATH, FNAMES, CENTER_LON, CENTER_LAT, THETA )
%   Calculate a wind-aligned line density for a given time period.
%
%   Calculates a line density of NO2 up and downwind of a city by aligning
%   each day's plume to the x-axis as described in Valin 2013.  By fitting
%   an exponentially modified Gaussian function to the line density,
%   certain features of the NOx emissions and chemistry can be derived.
%
%   c.f.    Beirle et al., Science, 2015, pp. 1737-1739
%           Valin et al., Geophys. Res. Lett., 2013, pp. 1856-1860
%           de Foy et al., Atmos. Environ., 2014, pp. 66-77
%           Lu et al., ACP, 2015, pp. 10367-10383
%
%   Required inputs:
%
%       fpath - the path to where the files containing a Data structure are
%       located. This is the structure in files output from BEHR_main.m and
%       must be acceptable as the first input to rotate_plume.m
%
%       fnames - any one of several methods of specifying the files to
%       read. Could be a structure output from the dir() command, a
%       string that when passed as dir(fullfile(fpath,fnames)) will return
%       the file(s) desired, or a cell array of file names.
%
%       center_lon, center_lat - the coordinates considered the center of
%       the domain, usually the coordinate of a city or other point
%       emission source.
%
%       theta - a vector of wind directions given as degrees CCW from east
%       between -180 and 180. Must match the number of files to be loaded.
%
%   Outputs:
%
%       no2_x - the x-coordinates of the line density (in km from center
%       lon/lat)
%
%       no2_linedensity - the line density in mol/km.
%
%       no2_linedens_std - the standard deviation of the line density,
%       calculated by first calculated a standard deviation of the column
%       densities across the time average weighted by the grid cell
%       areaweight, then adding up those in quadrature along the line of
%       integration, multiplied by the width of the grid cells.
%
%       lon, lat - the longitude and latitude coordinates of the
%       wind-aligned column densities.
%
%       no2_mean - the average wind-aligned column densities.
%
%       no2_std - the area-weighted standard deviation of the wind aligned
%       column densities.
%
%       num_valid_obs - an array the same size as no2_mean that counts the
%       number of observations that went into each column density. NaNs and
%       cases where areaweight = 0 do not count.
%
%       nox - 3D matrix of individual rotated swaths, mainly for debugging
%       purposes.
%
%       debug_cell - a cell array with the file name and swath number
%       corresponding to each 2D slice of the output NOX.
%
%   Parameter inputs:
%
%
%       data_ind - which top-level indices of the Data structure loaded from
%       the files to use. Can be a single index, a vector of indicies, or
%       (as is default) an empty vector which will use all the swaths present
%       in each file.
%
%       'windvel' - a vector of wind velocities to be used in filtering out
%       days that do not meet desired criteria. If not given, all days that
%       have sufficient coverage (pixels not removed for clouds or row
%       anomaly) are used.
%
%       'windcrit' - the number to compare windvel values to, if given, the
%       parameter 'windop' must also be specified.
%
%       'windop' - can be '<', '>', '<=', or '>=' and will be used with
%       'windcrit' to evalute windvel values in the expression windvel
%       <windop> windcrit. If that evaluates to false, the day will be
%       skipped.
%
%       'crit_logical' - an alternative to the wind inputs, this should be a 
%       logical vector that is true for days that shold be included. Must
%       have the same number of elements as the number of files to use.
%
%       'rel_box_corners' - a four element vector to be passed to rotate
%       plume describing how large a box to use to circumscribe the plumes.
%       See rotate_plume for more information.
%
%       'no_reject' - boolean, defaults to false. If true, skips rejecting
%       pixels for cloud fraction, row anomaly, etc. Only recommended if
%       the input data is model data.
%
%       'force_calc' - if true will override the criteria that rejects days
%       with too many unfilled pixels and use all days that meet the
%       windvel criterion. Defaults to false.
%
%       'days_of_week' - string understood as the DAYS_OF_WEEK argument to
%       DO_KEEP_DAY_OF_WEEK() specifying which days to keep. Default is
%       UMTWRFS, i.e. all days.
%
%       'wind_dir_weights' - a vector that gives a factor that VCDs from a
%       certain wind direction should be weighted by in the average.
%       Requires 'wind_weights_bins' as well. If not given, defaults to 1
%       for all directions.
%
%       'wind_weights_bins' - A vector giving the edges of the bins for the
%       weights given by 'wind_dir_weights'. If wind_dir_weights has N
%       elements, this must have N+1. Must be in the same units as THETA.
%
%       'linedens_field' - what fields in the Data structure loaded to
%       grid. Default is 'behr', does standard BEHR gridding. For other
%       options, see ROTATE_PLUME's "grid_control" parameter.
%
%       'grid_method' - see ROTATE_PLUME.
%
%       'DEBUG_LEVEL' - level of output to console. Defaults to 1, 0 shuts
%       off all output.
%
%   Josh Laughner <joshlaugh5@gmail.com> 5 Feb 2016

%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%% INPUT PARSING %%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%

E = JLLErrors;

p=inputParser;
p.addOptional('nox_or_no2','no2',@(x) ismember(lower(x),{'nox','no2'}));
p.addParameter('rel_box_corners',[]);
p.addParameter('no_reject',false);
p.addParameter('force_calc',false);
p.addParameter('days_of_week', 'UMTWRFS');
p.addParameter('wind_dir_weights',[]);
p.addParameter('wind_weights_bins',[]);
p.addParameter('linedens_field', 'behr');
p.addParameter('grid_method', 'cvm');
p.addParameter('DEBUG_LEVEL',1);

p.parse(varargin{:});

pout=p.Results;
nox_or_no2 = pout.nox_or_no2;
rel_box_corners = pout.rel_box_corners;
force_calc = pout.force_calc;
no_reject = pout.no_reject;
days_of_week = pout.days_of_week;
wind_dir_weights = pout.wind_dir_weights;
wind_weights_bins = pout.wind_weights_bins;
linedens_field = pout.linedens_field;
grid_method = pout.grid_method;
DEBUG_LEVEL = pout.DEBUG_LEVEL;

if ~ischar(fpath)
    E.badinput('fpath must be a string')
elseif ~exist(fpath,'dir')
    E.badinput('fpath (%s) does not exist',fpath);
end

if ischar(fnames)
    fnames_struct = dir(fullfile(fpath,fnames));
elseif iscellstr(fnames)
    fnames_struct = struct('name', repmat({''},size(fnames)));
    for f=1:numel(fnames)
        fnames_struct(f).name = fnames{f};
    end
elseif isstruct(fnames) && isfield(fnames,'name')
    fnames_struct = fnames;
else
    E.badinput('Input fnames is not a valid format; see documentation')
end

if ~isscalar(center_lon) || ~isnumeric(center_lon)
    E.badinput('center_lon must be a scalar number')
end
if ~isscalar(center_lat) || ~isnumeric(center_lat)
    E.badinput('center_lat must be a scalar number')
end

if ~isnumeric(theta) || any(theta(:) < -180 | theta(:) > 180) || size(theta,1) ~= numel(fnames_struct)
    E.badinput('theta must be a numeric vector with values between -180 and +180 that has the same number of elements as the number of files to be loaded')
end

if ~isequal(size(wind_logical), size(theta)) || ~islogical(wind_logical)
    E.badinput('wind_logical must be a logical matrix with the same size as THETA');
end

if ~isempty(rel_box_corners) && ( ~isvector(rel_box_corners) || numel(rel_box_corners) ~= 4 )
    E.badinput('rel_box_corners must be a 4 element vector');
end
if ~isscalar(force_calc) || ~islogical(force_calc)
    E.badinput('force_calc must be a scalar logical');
end
if ~isscalar(DEBUG_LEVEL) || ~isnumeric(DEBUG_LEVEL) || DEBUG_LEVEL < 0
    E.badinput('DEBUG_LEVEL must be a scalar number >= 0.')
end

if xor(isempty(wind_dir_weights), isempty(wind_weights_bins))
    E.badinput('Both or neither of ''wind_dir_weights'' and ''wind_weights_bins'' must be given')
elseif isempty(wind_dir_weights) % if we passed the xor() test if one is empty, both are
    wind_dir_weights = 1;
    wind_weights_bins = [-Inf Inf];
elseif ~isnumeric(wind_dir_weights) || ~isnumeric(wind_weights_bins)
    E.badinput('''wind_dir_weights'' and ''wind_weights_bins'' must both be numeric arrays')
elseif numel(wind_dir_weights) ~= numel(wind_weights_bins) - 1
    E.badinput('''wind_weights_bins'' must have one more element than ''wind_dir_weights''');
elseif any(diff(wind_weights_bins(:)) < 0)
    E.badinput('''wind_weights_bins'' must be monotonically increasing');
end

if ischar(linedens_field) && ~strcmpi(linedens_field, 'behr')
    % rotate_plume, if given 'behr' as grid_control will do the standard
    % rotation and gridding for BEHR VCDs. Otherwise it expects a cell
    % array of field names to grid.
    omi_field = linedens_field;
    linedens_field = {linedens_field};
elseif ~ischar(linedens_field)
    E.badinput('"linedens_field" must be a char array');
else
    omi_field = 'BEHRColumnAmountNO2Trop';
end

% Pixel reject structure for row anomaly, cloud fraction, etc.
reject_details = struct('cloud_type', 'omi', 'cloud_frac', 0.2, 'row_anom_mode', 'XTrackFlags', 'check_behr_amf', true);

%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%% MAIN FUNCTION %%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%

% Constants
if strcmpi(nox_or_no2,'nox')
    nox_no2_scale = 1.32; % scales NO2 to NOx columns, c.f. supporting info for Beirle et al. 2011 (Science)
else
    nox_no2_scale = 1;
end

%fid = fopen('index_date_and_swath.txt','w');

create_array = true;
i = 0;
for d=1:numel(fnames_struct)
    this_file = fullfile(fpath,fnames_struct(d).name);
    % Some days are not produced in BEHR. Those days need to be skipped.
    if exist(this_file, 'file')
        fprintf('Loading %s\n', this_file);
        D = load(this_file,'Data');
    else
        fprintf('%s does not exist, skipping\n', this_file);
        continue
    end
    
    if ~strcmpi(linedens_field, 'behr') && ~isfield(D.Data, linedens_field)
        E.callError('data_missing_field', 'File %s is missing the requested data field "%s"', this_file, linedens_field);
    end

    if ~do_keep_day_of_week(D.Data(1).Date, days_of_week)
        if DEBUG_LEVEL > 0
            fprintf('%s, day of week is not in %s, skipping\n', D.Data(1).Date, days_of_week);
        end
        continue
    end
    
    n_swath = numel(D.Data);
    for s=1:n_swath
        %fprintf(fid,'** %d: %s swath %d\n', i, fnames_struct(d).name, s);
        if ~wind_logical(d,s)
            if DEBUG_LEVEL > 0
                fprintf('Condition criteria not met, skipping %s orbit %d\n', fnames_struct(d).name, s);
            end
            continue
        end
        
        if isnan(theta(d,s))
            E.notimplemented('No handling for NaN wind direction')
        end
        
        if DEBUG_LEVEL > 0; disp('Rotating plume'); end
        if ~isempty(rel_box_corners)
            OMI = rotate_plume(D.Data(s), center_lon, center_lat, theta(d,s), rel_box_corners, 'grid_control', linedens_field, 'grid_method', grid_method);
        else
            OMI = rotate_plume(D.Data(s), center_lon, center_lat, theta(d,s), 'grid_control', linedens_field, 'grid_method', grid_method);
        end
        if isempty(OMI.Longitude)
            if DEBUG_LEVEL > 0; fprintf('No grid cells in %s\n',fnames_struct(d).name); end 
            continue
        end
        
        i = i+1;
        
        if no_reject
            if DEBUG_LEVEL > 1
                fprintf('  Not rejecting pixels\n');
            end
            xx = true(size(OMI.Longitude)); 
        else
            if DEBUG_LEVEL > 1
                fprintf('  Rejecting pixels\n');
            end
            OMI = omi_pixel_reject(OMI, 'detailed', reject_details);
            xx = OMI.Areaweight > 0;
        end
        
        if create_array
            nox = nan(size(OMI.Longitude,1), size(OMI.Longitude,2), numel(fnames_struct)*n_swath);
            aw = nan(size(OMI.Longitude,1), size(OMI.Longitude,2), numel(fnames_struct)*n_swath);
            lon = OMI.Longitude;
            lat = OMI.Latitude;
            debug_cell = cell(numel(fnames_struct)*n_swath,1);
            create_array = false;
        end
        
        debug_cell{i} = sprintf('%s: swath %d', fnames_struct(d).name, s);

            
        OMI.(omi_field)(~xx) = nan;
        OMI.Areaweight(~xx) = nan;
        % Multiply the weights by the wind direction weight, so that it
        % weights the NOx average and gets normalized out along with
        % the areaweight. The idea is that since I want to use slow
        % wind speeds to generate a source function as in Liu et al.
        % 2016 (doi: 10.5194/acp-16-5283-2016) but do it with rotated
        % line densities, I should weight the slow line densities so
        % that the different wind directions contribute to the line
        % density in the same proportions as in the fast line
        % densities.
        nox(:,:,i) = OMI.(omi_field) * nox_no2_scale;
        aw(:,:,i) = OMI.Areaweight * get_wind_dir_weight(theta(d,s), wind_dir_weights, wind_weights_bins);
        
    end
end
%fclose(fid);

%no2_mean = nanmean(nox,3);
no2_mean = nansum2(nox .* aw, 3) ./ nansum2(aw, 3);
num_valid_obs = sum( ~isnan(nox) & aw > 0, 3);
% Calculate the weighted standard deviation (c.f.
% https://en.wikipedia.org/wiki/Mean_square_weighted_deviation)
no2_var = (nansum2(aw .* nox.^2, 3) .* nansum2(aw,3) - (nansum2(aw .* nox, 3)).^2)./(nansum2(aw,3).^2 - nansum(aw.^2,3));
no2_std = sqrt(no2_var);

% Calculate the line density. See de Foy et al., Atmos. Environ. (2014) p.
% 66. Basically an integration along the line perpendicular to the plume.

if DEBUG_LEVEL > 0; disp('Calculating line density'); end
no2_linedens = zeros(1,size(lon,2));
no2_lindens_std = zeros(1,size(lon,2));
no2_x = nan(1,size(lon,2));
d_cm = nan(size(lon));
for a=1:numel(no2_linedens)
    for b=1:size(lon,1)-1 % because we need a latitudinal difference to average over, we can't do the last row.
        d_cm(b,a) = m_lldist(lon(b:b+1,a),lat(b:b+1,a)) * 1e5;
        if ~isnan(no2_mean(b,a))
            no2_linedens(a) = no2_linedens(a) + no2_mean(b,a) * d_cm(b,a);
            % add the uncertainties in quadrature.
            no2_lindens_std(a) = no2_lindens_std(a) + no2_std(b,a).^2 * d_cm(b,a);
        end
    end
    % Calculate x in km distant from the center lon/lat. OMI is a gridded
    % representation so OMI.Longitude(:,a) are all the same, as are
    % OMI.Latitude(b,:).
    no2_x(a) = m_lldist([lon(1,a), center_lon], [center_lat, center_lat]) * sign(lon(1,a) - center_lon);
end

% Remove any values that never got anything added to them b/c there were no
% non-nan values for that transect.
rr = ~all(isnan(no2_mean),1);
no2_x = no2_x(rr);
no2_linedens = no2_linedens(rr);
no2_lindens_std = no2_lindens_std(rr);

% Finalize the uncertainties
no2_lindens_std = sqrt(no2_lindens_std);

% Convert line density from molec/cm to moles/km, if doing standard BEHR
% line densities
if strcmpi(linedens_field, 'behr')
    no2_linedens = no2_linedens * 1e5 / 6.022e23;
    no2_lindens_std = no2_lindens_std * 1e5 / 6.022e23;
end
end


function weight = get_wind_dir_weight(this_theta, wind_dir_weights, wind_dir_weight_bins)
idx = find(this_theta > wind_dir_weight_bins, 1, 'last');
if isempty(idx)
    E = JLLErrors;
    E.callError('undefined_wind_dir_weight', 'No wind direction weight defined for theta = %.1f', this_theta);
end
weight = wind_dir_weights(idx);
end
