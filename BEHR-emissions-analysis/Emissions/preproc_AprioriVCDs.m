function [] = preproc_AprioriVCDs(start_date, end_date, save_path, varargin)
%UNTITLED2 Summary of this function goes here
%   Detailed explanation goes here

E = JLLErrors;

p = inputParser;
p.addParameter('region','us');
p.addParameter('prof_mode','daily');
%p.addParameter('behr_mat_dir','');
p.addParameter('overwrite',false);
p.addParameter('n_workers', 0);
p.addParameter('DEBUG_LEVEL',2);

p.parse(varargin{:});
pout = p.Results;

start_date = validate_date(start_date);
end_date = validate_date(end_date);

region = pout.region;
prof_mode = pout.prof_mode;
%behr_mat_dir = pout.behr_mat_dir;
do_overwrite = pout.overwrite;
n_workers = pout.n_workers;
DEBUG_LEVEL = pout.DEBUG_LEVEL;

if ~ischar(save_path)
    E.badinput('SAVE_PATH must be a char array');
elseif ~exist(save_path, 'dir')
    E.badinput('SAVE_PATH (%s) is not a valid directory', save_path);
end

% if ~ischar(behr_mat_dir)
%     E.badinput('The parameter "behr_mat_dir" must be a char array');
% elseif ~exist(behr_mat_dir, 'dir')
%     E.badinput('The path for parameter "behr_mat_dir" (%s) does not exist', behr_mat_dir);
% elseif isempty(behr_mat_dir)
%     behr_mat_dir = behr_paths.BEHRMatSubdir(region, prof_mode);
% end

if ~islogical(do_overwrite) || ~isscalar(do_overwrite)
    E.badinput('The parameter "overwrite" must be a scalar logical value');
end

%%%%%%%%%%%%%%%%%
% MAIN FUNCTION %
%%%%%%%%%%%%%%%%%

dvec = make_dvec(start_date, end_date);
req_fields = get_required_fields();
behr_mat_dir = behr_paths.BEHRMatSubdir(region, prof_mode);

parfor (i_day = 1:numel(dvec), n_workers)
    if DEBUG_LEVEL > 0
        fprintf('Working on %s\n', datestr(dvec(i_day)));
    end
    full_save_name = fullfile(save_path, make_save_name(dvec(i_day), prof_mode, region));
    if exist(full_save_name, 'file') && ~do_overwrite
        fprintf('Output file %s already exists, skipping\n', full_save_name);
        continue
    end

    file_to_load = fullfile(behr_mat_dir, behr_filename(dvec(i_day), prof_mode, region));
    if ~exist(file_to_load, 'file')
        fprintf('%s does not exist, skipping\n', file_to_load);
        continue
    end

    if DEBUG_LEVEL > 1
        fprintf('  Loading %s\n', file_to_load);
    end
    LoadTmp = load(file_to_load,'Data');
    DataBEHR = LoadTmp.Data;
    
    Data = make_empty_struct_from_cell(req_fields);
    Data = repmat(Data, size(DataBEHR));
    
    for i_orbit = 1:numel(Data)
        if DEBUG_LEVEL > 1
            fprintf('    Orbit %d of %d\n', i_orbit, numel(Data));
        end
        Data(i_orbit) = copy_structure_fields(DataBEHR(i_orbit),Data(i_orbit));
        Data(i_orbit).BEHRColumnAmountNO2Trop = nan(size(Data(i_orbit).BEHRColumnAmountNO2Trop));
        for i_pix = 1:numel(Data(i_orbit).Longitude)
            Data(i_orbit).BEHRColumnAmountNO2Trop(i_pix) = integrate_wrf_profile(DataBEHR(i_orbit).BEHRNO2apriori(:,i_pix),...
                DataBEHR(i_orbit).BEHRPressureLevels(:,i_pix), DataBEHR(i_orbit).BEHRSurfacePressure(i_pix), DataBEHR(i_orbit).BEHRTropopausePressure(i_pix));
        end
    end
    
    if DEBUG_LEVEL > 0
        fprintf('  Saving as %s\n', full_save_name);
    end
    save_helper(full_save_name, Data);
end

end

function save_helper(save_name, Data)
save(save_name, 'Data');
end

function vcd = integrate_wrf_profile(no2_prof, pres_levs, surf_pres, trop_pres)
xx = ~isnan(pres_levs);
if sum(xx) < 2
    vcd = nan;
    return;
end

if surf_pres > max(pres_levs(xx)) && surf_pres <= max(behr_pres_levels())
    error('wrf_prof:nan_pres', 'Surface pressure is below the first non-NaN profile level');
end

vcd = integPr2(no2_prof(xx), pres_levs(xx), surf_pres, trop_pres);

end

function dvec = make_dvec(start_dates, end_dates)
dvec = [];
for a=1:numel(start_dates)
    dvec = veccat(dvec, start_dates(a):end_dates(a));
end
end

function fields = get_required_fields()
behr_paths.SetPythonPath();
fields = python2matlab(py.PSM_Main.behr_datasets('cvm'));
% Date is needed by calc_line_density to pick day of week;
% ViewingZenithAngle is needed by rotate_plume to remove grid cells at too
% extreme a viewing angle. Swath is needed in psm_wrapper.
fields = veccat(fields, {'Date', 'ViewingZenithAngle','Swath'},'column');
end

function save_name = make_save_name(curr_date, prof_mode, region)
save_name = behr_filename(curr_date, prof_mode, region);
save_name = strrep(save_name, 'OMI', 'WRF');
end
