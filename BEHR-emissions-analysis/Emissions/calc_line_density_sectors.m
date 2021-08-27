function [ no2_x, no2_linedens, no2_lindens_std, lon, lat, no2_mean, no2_std, num_valid_obs, nox, debug_cell ] = calc_line_density_sectors( fpath, fnames, center_lon, center_lat, theta, wind_logical, varargin )
%[ NO2_X, NO2_LINEDENS, NO2_LINEDENS_STD, LON, LAT, NO2_MEAN, NO2_STD, NUM_VALID_OBS] = CALC_LINE_DENSITY( FPATH, FNAMES, CENTER_LON, CENTER_LAT, THETA )
%   Calculate a wind-aligned line density for a given time period.
%
%   Calculates a line density of NO2 up and downwind of a city for 8
%   different wind sectors a la Beirle et al. 2011.
%
%   c.f.    Beirle et al., Science, 2015, pp. 1737-1739
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
%       theta - a matrix of wind directions given as degrees CCW from east
%       between -180 and 180. The first dimension must match the number of
%       files to be loaded; the second must correspond to the orbit index
%       in each Data structure; i.e. the winds for the third orbit,
%       Data(3), in the 10th file would be given in theta(10,3).
%
%       wind_logical - used to separate sufficiently fast winds, this
%       should be a logical matrix that is true for orbits that shold be
%       included. Must have the same configuration as THETA.
%
%   Outputs:
%
%       no2_x - the x-coordinates of the line density (in km from center
%       lon/lat) as a structure for each direction.
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
%       windvel criterion. Defaults to false; only used if 'interp' is
%       true.
%
%       'interp' - boolean, defaults to false (CHANGED 6 Jul 2018). If
%       true, individual days are chosen if they have a sufficient number
%       of viable observations to represent the entire domain well. Any
%       missing elements are filled my interpolation.  If false, bad
%       observations (cloud fraction or row anomaly) are removed, but all
%       days are averaged together. This mode is closer to how Lu et al.
%       2015 did their analysis, I believe.
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
p.addParameter('force_calc',false);
p.addParameter('interp',false);
p.addParameter('no_reject',false);
p.addParameter('DEBUG_LEVEL',1);

p.parse(varargin{:});

pout=p.Results;
nox_or_no2 = pout.nox_or_no2;
rel_box_corners = pout.rel_box_corners;
force_calc = pout.force_calc;
interp_bool = pout.interp;
no_reject = pout.no_reject;
DEBUG_LEVEL = pout.DEBUG_LEVEL;

if ~ischar(fpath)
    E.badinput('fpath must be a string')
elseif ~exist(fpath,'dir')
    E.badinput('fpath (%s) does not exist',fpath);
end

if ischar(fnames)
    fnames_struct = dir(fullfile(fpath,fnames));
elseif iscellstr(fnames);
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
    E.badinput('theta must be a numeric vector with values between -180 and +180 that has the same size in the first dimension as the number of files to be loaded')
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

if ~isscalar(interp_bool) || ~islogical(interp_bool)
    E.badinput('interp must be a scalar logical')
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

create_array = true;
i = 0;
for d=1:numel(fnames_struct)
    D = load(fullfile(fpath,fnames_struct(d).name),'Data');

    % We'll still "rotate" each day but divide it by sectors. This isn't so
    % much trying to reproduce Beirle 11 as give me a way to find out which
    % directions contribute to certain features of the line density.
    n_swath = numel(D.Data);
    for s=1:n_swath
        if ~wind_logical(d,s)
            if DEBUG_LEVEL > 0; 
                fprintf('Condition criteria not met, skipping %s orbit %d\n',fnames_struct(d).name, s); 
            end
            continue
        end
        
        
        if DEBUG_LEVEL > 0; disp('Rotating plume'); end
        if ~isempty(rel_box_corners)
            OMI = rotate_plume(D.Data(s), center_lon, center_lat, theta(d,s), rel_box_corners);
        else
            OMI = rotate_plume(D.Data(s), center_lon, center_lat, theta(d,s));
        end
        if isempty(OMI.Longitude)
            if DEBUG_LEVEL > 0; fprintf('No grid cells in %s\n',fnames_struct(d).name); end 
            continue
        end
        
        i = i+1;
        
        if no_reject
            xx = true(size(OMI.Longitude));
        else
            OMI = omi_pixel_reject(OMI, 'detailed', reject_details);
            xx = OMI.Areaweight > 0;
        end
        if create_array
            directions = {'W','SW','S','SE','E','NE','N','NW'};
            theta_bin_edges = [-180, -157.5, -112.5, -67.5, -22.5, 22.5, 67.5, 112.5, 157.5];
            nox = make_empty_struct_from_cell(directions, nan(size(OMI.Longitude,1), size(OMI.Longitude,2), numel(fnames_struct)*n_swath));
            aw = make_empty_struct_from_cell(directions, nan(size(OMI.Longitude,1), size(OMI.Longitude,2), numel(fnames_struct)*n_swath));
            lon = OMI.Longitude;
            lat = OMI.Latitude;
            debug_cell = cell(numel(fnames_struct)*n_swath,1);
            create_array = false;
        end
        
        debug_cell{i} = sprintf('%s: swath %d', fnames_struct(d).name, s);

        
        % Find which bin this belongs in
        if theta(d,s) < theta_bin_edges(2) || theta(d,s) > theta_bin_edges(end)
            bin = 'W';
        else
            theta_xx = theta_bin_edges(1:end-1) <= theta(d,s) & theta_bin_edges(2:end) > theta(d,s);
            bin = directions{theta_xx};
        end
        
        if interp_bool
            % This criterion accounts for how many neighbors are empty, giving more
            % weight to large clumps of NaNs (due to row anomaly or clouds) and
            % less weight to scattered pixels. It still looks like 50% is a good
            % cutoff.
            mfrac = badpix_metric(~xx);
            if mfrac > .5 && ~force_calc
                if DEBUG_LEVEL > 0; fprintf('Too many clumped missing pixels, skipping %s\n',fnames_struct(d).name); end
                continue
            end
            
            % Fill in empty grid boxes by interpolation to prevent discontinuity in the
            % line density due to an uneven number of missing values for different
            % distances from the city.
            if DEBUG_LEVEL > 0; disp('Interpolating to fill in gaps in NO2 matrix'); end
            F = scatteredInterpolant(OMI.Longitude(xx), OMI.Latitude(xx), OMI.BEHRColumnAmountNO2Trop(xx));
            Faw = scatteredInterpolant(OMI.Longitude(xx), OMI.Latitude(xx), OMI.Areaweight(xx));
            nox.(bin)(:,:,i) = F(OMI.Longitude, OMI.Latitude)*nox_no2_scale;
            aw.(bin)(:,:,i) = Faw(OMI.Longitude, OMI.Latitude);
        else
            
            OMI.BEHRColumnAmountNO2Trop(~xx) = nan;
            OMI.Areaweight(~xx) = nan;
            nox.(bin)(:,:,i) = OMI.BEHRColumnAmountNO2Trop*nox_no2_scale;
            aw.(bin)(:,:,i) = OMI.Areaweight;
        end
    end
end

%no2_mean = nanmean(nox,3);
for f=1:numel(directions)
    no2_mean.(directions{f}) = nansum2(nox.(directions{f}) .* aw.(directions{f}), 3) ./ nansum2(aw.(directions{f}), 3);
    num_valid_obs = sum( ~isnan(nox.(directions{f})) & aw.(directions{f}) > 0, 3);
    % Calculate the weighted standard deviation (c.f.
    % https://en.wikipedia.org/wiki/Mean_square_weighted_deviation)
    no2_var.(directions{f}) = (nansum2(aw.(directions{f}) .* nox.(directions{f}).^2, 3) .* nansum2(aw.(directions{f}),3) - (nansum2(aw.(directions{f}) .* nox.(directions{f}), 3)).^2)./(nansum2(aw.(directions{f}),3).^2 - nansum(aw.(directions{f}).^2,3));
    no2_std.(directions{f}) = sqrt(no2_var.(directions{f}));
    
    % Calculate the line density. See de Foy et al., Atmos. Environ. (2014) p.
    % 66. Basically an integration along the line perpendicular to the plume.
    
    if DEBUG_LEVEL > 0; disp('Calculating line density'); end
    no2_linedens.(directions{f}) = zeros(1,size(lon,2));
    no2_lindens_std.(directions{f}) = zeros(1,size(lon,2));
    no2_x.(directions{f}) = nan(1,size(lon,2));
    d_cm = nan(size(lon));
    for a=1:numel(no2_linedens.(directions{f}))
        for b=1:size(lon,1)-1 % because we need a latitudinal difference to average over, we can't do the last row.
            d_cm(b,a) = m_lldist(lon(b:b+1,a),lat(b:b+1,a)) * 1e5;
            if ~isnan(no2_mean.(directions{f})(b,a))
                no2_linedens.(directions{f})(a) = no2_linedens.(directions{f})(a) + no2_mean.(directions{f})(b,a) * d_cm(b,a);
                % add the uncertainties in quadrature.
                no2_lindens_std.(directions{f})(a) = no2_lindens_std.(directions{f})(a) + no2_std.(directions{f})(b,a).^2 * d_cm(b,a);
            end
        end
        % Calculate x in km distant from the center lon/lat. OMI is a gridded
        % representation so OMI.Longitude(:,a) are all the same, as are
        % OMI.Latitude(b,:).
        no2_x.(directions{f})(a) = m_lldist([lon(1,a), center_lon], [center_lat, center_lat]) * sign(lon(1,a) - center_lon);
    end
    
    % Remove any values that never got anything added to them b/c there were no
    % non-nan values for that transect.
    rr = ~all(isnan(no2_mean.(directions{f})),1);
    no2_x.(directions{f}) = no2_x.(directions{f})(rr);
    no2_linedens.(directions{f}) = no2_linedens.(directions{f})(rr);
    no2_lindens_std.(directions{f}) = no2_lindens_std.(directions{f})(rr);
    
    % Finalize the uncertainties
    no2_lindens_std.(directions{f}) = sqrt(no2_lindens_std.(directions{f}));
    
    % Convert line density from molec/cm to moles/km
    no2_linedens.(directions{f}) = no2_linedens.(directions{f}) * 1e5 / 6.022e23;
    no2_lindens_std.(directions{f}) = no2_lindens_std.(directions{f}) * 1e5 / 6.022e23;
end
end

function [mfrac, msum] = badpix_metric(xx)
% Calculates a metric for the badness of missing pixels. A missing grid
% cell doesn't contribute to this metric unless it has at least one
% neighbor that is also missing
sz = size(xx);
m = zeros(sz);
for a=1:sz(1)
    for b=1:sz(2)
        if a > 1 && b > 1
            m(a,b) = m(a,b) + xx(a-1, b-1);
        end
        if a > 1
            m(a,b) = m(a,b) + xx(a-1,b);
        end
        if a > 1 && b < sz(2)
            m(a,b) = m(a,b) + xx(a-1,b+1);
        end
        if b < sz(2)
            m(a,b) = m(a,b) + xx(a,b+1);
        end
        if a < sz(1) && b < sz(2)
            m(a,b) = m(a,b) + xx(a+1,b+1);
        end
        if a < sz(1)
            m(a,b) = m(a,b) + xx(a+1,b);
        end
        if a < sz(1) && b > 1
            m(a,b) = m(a,b) + xx(a+1,b-1);
        end
        if b > 1
            m(a,b) = m(a,b) + xx(a,b-1);
        end
    end
end

% How many neighbors there are, i.e. the maximum value m could take on.
% Corner cells have 3 neighbors, non-corner edge cells have 5 and non-edge
% cells have 8.
n_neighbors = prod(sz-2)*8 + (prod(sz) - prod(sz-2) - 4) * 5 + 12;
msum = sum(m(:));
mfrac = msum / n_neighbors;
end
