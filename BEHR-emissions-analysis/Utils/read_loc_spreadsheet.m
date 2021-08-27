function [ locs ] = read_loc_spreadsheet(  )
%UNTITLED2 Summary of this function goes here
%   Detailed explanation goes here

mydir = fileparts(mfilename('fullpath'));
spreadsheet_file = fullfile(mydir, '..', 'Workspaces', 'SiteData', 'trend_locations.xlsx');
[~,sheets] = xlsfinfo(spreadsheet_file);

locs = [];

for a=1:numel(sheets)
    [~,~,raw] = xlsread(spreadsheet_file, sheets{a});
    % Only take columns with headers, some rows have unneeded notes at the
    % end
    xx = cellfun(@(x) ischar(x), raw(1,:)); 
    
    % Stop at a blank line - some sheets have notes at the end
    nans = cellfun(@(x) ~ischar(x) && isnan(x), raw);
    last_line = find(all(nans,2),1) - 1;
    if isempty(last_line)
        last_line = size(raw,1);
    end
    
    header = raw(1,xx);
    raw = raw(2:last_line,xx);
    
    % Convert the box size from a string to a numeric array
    box_size_ind = strcmpi(header, 'BoxSize');
    if sum(box_size_ind) == 0
        error('trend_locations:missing_category', 'The category "BoxSize" is not present in the Excel file')
    else
        raw(:,box_size_ind) = convert_box_size(raw(:,box_size_ind), sheets{a});
    end
    
    % Parse the WindRejects field so that it is an N-by-2 array giving each
    % range to reject as a list of N ranges of wind headings, 
    wind_reject_ind = strcmpi(header, 'WindRejects');
    if sum(wind_reject_ind) == 0
        error('trend_locations:missing_category', 'The category "WindRejects" is not present in the Excel file')
    else
        raw(:,wind_reject_ind) = parse_wind_rejects(raw(:,wind_reject_ind), sheets{a});
    end
    
    % Parse the WindRejects field so that it is an N-by-2 array giving each
    % range to reject as a list of N ranges of wind headings, 
    wrf_wind_reject_ind = strcmpi(header, 'WRFWindRejects');
    if sum(wrf_wind_reject_ind) == 0
        error('trend_locations:missing_category', 'The category "WRFWindRejects" is not present in the Excel file')
    else
        raw(:,wrf_wind_reject_ind) = parse_wind_rejects(raw(:,wrf_wind_reject_ind), sheets{a});
    end
    
    
    header{end+1} = 'SiteType';
    raw(:,end+1) = repmat({sheets{a}},size(raw,1),1);
    
    these_locs = cell2struct(raw,header,2);
    tmp_locs = cat(1, locs, these_locs);
    locs = tmp_locs;
end

end

function cell_out = convert_box_size(cell_in, location_type)
cell_out = cell(size(cell_in));
for a=1:numel(cell_in)
    if isnan(cell_in{a})
        if strcmpi(location_type, 'cities')
            cell_out{a} = [1 2 1 1];
        elseif strcmpi(location_type, 'powerplants')
            cell_out{a} = [0.5 1 0.5 0.5];
        elseif strcmpi(location_type, 'ruralareas')
            cell_out{a} = nan(1,4);
        else
            error('trend_locations:unknown_site_type','No default box size defined for site type "%s"', location_type);
        end
    else
        value = regexprep(cell_in{a},'[^\d\s\.]', '');
        cell_out{a} = str2double(strsplit(value));
    end
end
end

function reject_ranges_cell = parse_wind_rejects(cell_in, current_sheet)
sector.N = [67.5 112.5];
sector.NE = [22.5 67.5];
sector.E = [-22.5 22.5];
sector.SE = [-67.2 -22.5];
sector.S = [-112.5 -67.5];
sector.SW = [-152.5 -112.5];
sector.W = [152.5 -152.5];
sector.NW = [112.5 152.5];

sector_fns = fieldnames(sector);

reject_ranges_cell = cell(size(cell_in));
for a=1:numel(cell_in)
    if isnan(cell_in{a})
        reject_ranges_cell{a} = [];
    else
        % Split each cell into individual ranges, assuming the are comma
        % separated
        site_range_strings = split_on_commas_outside_brackets(cell_in{a});
        site_range_strings = cellfun(@strtrim, site_range_strings, 'UniformOutput', false);
        
        % For each range, see if it is one of the predefined sectors. If
        % not, then it must be a range of the form [ min, max ].
        site_ranges = nan(numel(site_range_strings), 2);
        
        for b=1:numel(site_range_strings)
            if ismember(site_range_strings{b}, sector_fns)
                site_ranges(b,:) = sector.(site_range_strings{b});
            else
                tmp_range = str2num(site_range_strings{b}); %#ok<ST2NM> not operating on scalar values - expecting two element vectors
                if isempty(tmp_range)
                    % Give the line number as a+1 to account for the header
                    error('trend_locations:wind_reject_parse_error', 'The wind reject value in line %d of %s has at least one value that cannot be interpreted as a sector or range', a+1, current_sheet);
                elseif numel(tmp_range) ~= 2
                    error('trend_locations:wind_reject_parse_error', 'One of the ranges given for wind reject in line %d of %s does not have exactly two elements', a+1, current_sheet);
                end
                
                site_ranges(b,:) = tmp_range;
            end
        end
        
        reject_ranges_cell{a} = site_ranges;
    end
end
end

function cell_out = split_on_commas_outside_brackets(str_in)
bracket_count = 0;
cell_out = {};
last_split_index = 1;
for a=1:length(str_in)
    if strcmp(str_in(a), ',') && bracket_count == 0
        cell_out{end+1} = str_in(last_split_index:a-1);
        last_split_index = a+1;
    elseif strcmp(str_in(a),'[')
        bracket_count = bracket_count + 1;
    elseif strcmp(str_in(a),']')
        bracket_count = bracket_count - 1;
    end
    
    if bracket_count < 0
        error('trend_locations:unmatched_bracket', 'Unmatched closing bracket found in string "%s"', str_in);
    end
end

cell_out{end+1} = str_in(last_split_index:end);
end