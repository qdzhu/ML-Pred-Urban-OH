classdef misc_prepare_ml_input
    
    methods(Static = true)
        
        function make_location_small_radius_wrf_avgs_file(time_periods, varargin)
            E = JLLErrors;
            p = advInputParser;
            p.addParameter('overwrite_vars', {});
            p.addParameter('read_qhcho',true);
            p.addParameter('skip_wrf', false);
            p.addParameter('loc_inds', 1:49);
            p.parse(varargin{:});
            pout = p.Results;
            overwrite_vars = pout.overwrite_vars;
            loc_inds = pout.loc_inds;
            read_qhcho = pout.read_qhcho;
            skip_wrf = pout.skip_wrf;
            %'hcho_vcd_aks',, 'no2_aks', 'hcho_aks'
            wrf_3d_vars = {'no', 'no2', 'ho', 'pres', 'hcho', 'iso','hcho_aks','qhcho_aks'};
            wrf_2d_vars = {'ndens', 'temperature', 'ho', 'LNOXHNO3', 'LNOXA','pan','no','no2','LNOXPAN','iso','o3',...
                'hcho','PHOTR_O31D','PHOTR_CH2OR','QVAPOR','PO3','pressure','COSZEN','no2_vcd','iso_vcd','hcho_vcd','PHOTR_NO2',...
                'PBLH','IC_FLASHCOUNT','CG_FLASHCOUNT','elev_bot','elev_top','hcho_vcd_aks','qhcho_vcd_aks'};
            wrf_vars = veccat(wrf_3d_vars, wrf_2d_vars);

            
            base_locs = misc_ml_oh_prepare.match_cbsa_to_cities();
            base_locs_names = {base_locs.ShortName};
            
            n_files_per_day = 6; % should match the assumed number for the winds files
            n_2d_vars = numel(wrf_2d_vars);
            n_3d_vars = numel(wrf_3d_vars);
            n_vars = numel(wrf_vars);
            n_wrf_levels = 29;
            n_locs = 49;
            n_times = numel(time_periods);
            
            winds_data = cell(n_times,1);
            exist_vars_flag = false(n_times, n_vars);
            
            for i_time = 1:n_times
                [start_date, end_date] = misc_emissions_analysis.select_start_end_dates(time_periods{i_time});
                if isfile(misc_emissions_analysis.wrf_small_radius_data_file_name(start_date, end_date))
                    winds_data{i_time} = load(misc_emissions_analysis.wrf_small_radius_data_file_name(start_date, end_date));
                    increment = true;
                else
                    winds_data{i_time} = load(misc_emissions_analysis.wrf_data_file_name(start_date, end_date));
                    increment = false;
                end
                if exist('loc_inds', 'var')
                    winds_data{i_time}.locs = misc_emissions_analysis.cutdown_locs_by_index(winds_data{i_time}.locs, loc_inds);
                end
                % Double check that the winds structs match the base locs
                tmp_locs = winds_data{i_time}.locs;
                loc_names = {tmp_locs.ShortName};
                
                % Append a substructure to put the WRF data into. Have one
                % value per day, we'll average together whichever hours are
                % used.
                field_2d = 'Averaged';
                field_3d = 'Profile';
                for i_loc=1:n_locs
        
                    for i_var = 1:numel(wrf_2d_vars)
                        data_struct_2d.(wrf_2d_vars{i_var}) = cell(1, size(tmp_locs(i_loc).WindUsedBool,1));
                    end
                    data_struct_2d.WindDir = cell(1, size(tmp_locs(i_loc).WindUsedBool,1));
                    data_struct_2d.WindSpeed = cell(1, size(tmp_locs(i_loc).WindUsedBool,1));
                    
                    for i_var = 1:numel(wrf_3d_vars)
                        data_struct_3d.(wrf_3d_vars{i_var}) = cell(1, size(tmp_locs(i_loc).WindUsedBool,1));
                    end
                    tmp_locs(i_loc).WRFData.(field_2d) = data_struct_2d;
                    tmp_locs(i_loc).WRFData.(field_3d) = data_struct_3d;
                    
                    for i_var = 1:n_vars
                        is_2d = is_var_2d(i_var);
                        var_name = wrf_vars{i_var};
                        if is_2d
                            substruct = field_2d;
                        else
                            substruct = field_3d;
                        end
                        if increment
                            exist_vars = fieldnames(winds_data{i_time}.locs(i_loc).WRFData.(substruct));

                            if any(strcmpi(exist_vars, var_name))
                                if any(strcmpi(overwrite_vars, var_name))
                                    fprintf('    Overwrite %s\n', wrf_vars{i_var});
                                else
                                    fprintf('    Reading %s from existing file\n', wrf_vars{i_var});
                                    tmp_locs(i_loc).WRFData.(substruct).(var_name)= winds_data{i_time}.locs(i_loc).WRFData.(substruct).(var_name);
                                    exist_vars_flag(i_time, i_var) = true;
                                end
                            end
                        end
                    end
                end
                
                
                winds_data{i_time}.locs = tmp_locs;
                
                winds_data{i_time}.savename = misc_emissions_analysis.wrf_small_radius_data_file_name(start_date, end_date);
            end
            
            all_years = unique(veccat(time_periods{:}));
            [total_starts, total_ends] = misc_emissions_analysis.select_start_end_dates(all_years);
            total_dvec = misc_emissions_analysis.make_datevec(total_starts, total_ends);
            
            WRF_Files_Getter = BEHRMatchedWRFFiles('DEBUG_LEVEL', 1);
            
            for d=1:numel(total_dvec)
                this_date = total_dvec(d);
                fprintf('Now on %s (%d of %d)\n', datestr(this_date), d, numel(total_dvec));
                % We need to pick the WRF file closest in time to the OMI
                % overpass, so we will load the BEHR file for this day and
                % calculate which WRF file is closest in time to that
                % swath, then calculate wind speed and direction for each
                % file. Later, we will actually pick the swath-specific
                % wind direction for each rotation.
                
                % If we're in the same month as the last time through this
                % loop, the closest_wrf method doesn't need to get the
                % directory listing of the WRF directory again (because it
                % should be organized by month and year).
                if skip_wrf
                    wrf_data_rd = load(misc_wrf_lifetime_analysis.wrf_avg_name(2014,'no2'));
                    wrf_lon = double(wrf_data_rd.xlon);
                    wrf_lat = double(wrf_data_rd.xlat);
                    fprintf('NO need to read wrf files \n');
                else
                    fprintf('%s: Gathering WRF files\n', datestr(total_dvec(d)));
                    wrf_files = WRF_Files_Getter.get_files_for_date(total_dvec(d));
                    if isempty(wrf_files)
                        continue
                    end

                    data_for_today = cell(n_wrf_levels, n_files_per_day, n_vars, n_locs);

                    for i_file=1:numel(wrf_files)
                        % Load the bottom five layers of each variable in turn,
                        % average it to the cities radii, then for each time
                        % period, figure out if this date is in it, if so, get
                        % the data for the right time for that city, and store
                        % it.
                        wrf_lon = ncread(wrf_files{i_file}, 'XLONG');
                        wrf_lat = ncread(wrf_files{i_file}, 'XLAT');
                        for i_var = 1:n_vars
                            if increment && exist_vars_flag(i_time, i_var)
                                %fprintf('The variable %s is already exist\n', wrf_vars{i_var});
                                continue;
                            end
                            if any(strfind(wrf_vars{i_var}, '_upwind')) && ~any(strfind(wrf_vars{i_var}, '_vcd'))
                                read_upwind = true;
                                fprintf('The variable %s is for upwind cells \n', wrf_vars{i_var});
                            else
                                read_upwind = false;
                            end

                            if any(strfind(wrf_vars{i_var}, '_surr')) && ~any(strfind(wrf_vars{i_var}, '_vcd'))
                                read_surr = true;
                                fprintf('The variable %s is for surrounding cells \n', wrf_vars{i_var});
                            else
                                read_surr = false;
                            end

                            if any(strfind(wrf_vars{i_var}, '_vcd')) || any(strfind(wrf_vars{i_var}, '_aks'))
                                fprintf('The variable %s is on hold \n', wrf_vars{i_var});
                                continue;
                            end

                            if strfind(wrf_vars{i_var}, 'FLASH')
                                read_flash = true;
                            else
                                read_flash = false;
                            end

                            if strfind(wrf_vars{i_var}, 'elev')
                                read_elev = true;
                            else
                                read_elev = false;
                            end

                            is_2d = is_var_2d(i_var);
                            if read_surr
                                var_name = strrep(wrf_vars{i_var}, '_surr', '');
                            elseif read_upwind
                                var_name = strrep(wrf_vars{i_var}, '_upwind', '');
                            elseif read_elev
                                var_name = 'elev';
                            else
                                var_name = wrf_vars{i_var};
                            end

                            fprintf('    Loading %s\n', wrf_vars{i_var});

                            if is_2d
                                num_lev = 5;
                            else
                                num_lev = Inf;
                            end

                            try
                                wrf_value = read_wrf_preproc(wrf_files{i_file}, var_name, [1 1 1 1], [Inf Inf num_lev Inf]);
                            catch err
                                if any(strcmpi(err.identifier, {'MATLAB:imagesci:netcdf:unableToOpenFileforRead','MATLAB:imagesci:netcdf:unknownLocation'}))
                                    fprintf('Cannot read %s from %s, try calculation by ourselves', wrf_vars{i_var}, wrf_files{i_file});
                                    try
                                        wrf_value = misc_wrf_lifetime_analysis.read_wrf_preproc_chem(wrf_files{i_file}, var_name, [1 1 1 1], [Inf Inf num_lev Inf]);
                                    catch
                                        if any(strcmpi(err.identifier, {'MATLAB:imagesci:netcdf:unableToOpenFileforRead','MATLAB:imagesci:netcdf:unknownLocation'}))
                                            fprintf('Fail to calculate %s from %s it will stay a NaN\n', wrf_vars{i_var}, wrf_files{i_file});
                                            continue
                                        else
                                            rethrow(err);
                                        end
                                    end
                                else
                                    if any(strcmpi(err.identifier, {'MATLAB:imagesci:netcdf:indexElementLength'}))
                                        wrf_value = read_wrf_preproc(wrf_files{i_file}, var_name, [1 1 1], [Inf Inf Inf]);
                                    else  
                                        rethrow(err)
                                    end
                                end 
                            end

                            if is_2d
                                if read_elev
                                    if strcmpi(wrf_vars{i_var}, 'elev_bot')
                                        wrf_value = wrf_value(:,:,1);
                                    else
                                        wrf_value = wrf_value(:,:,5);
                                    end
                                else
                                    wrf_value = nanmean(wrf_value, 3);
                                end
                            else
                                wrf_value = permute(wrf_value, [3 1 2]);
                                if size(wrf_value,1) ~= n_wrf_levels
                                    E.callError('wrf_level_mismatch', 'WRF levels for variable "%s" not the expected %d', var_name, n_wrf_levels);
                                end
                            end

                            if read_flash
                                this_hour = hour(datestr(date_from_wrf_filenames(wrf_files{i_file})));
                                prev_wrf_file = strrep(wrf_files{i_file}, sprintf('%d-00-00', this_hour), sprintf('%d-00-00', this_hour-1));
                                prev_value = read_wrf_preproc(prev_wrf_file, var_name, [1 1 1], [Inf Inf Inf]);
                                wrf_value = wrf_value - nanmean(prev_value,3);
                            end

                            for i_loc = 1:n_locs
                                if read_surr
                                    xx = misc_emissions_analysis.find_indices_in_radius_surr_loc(base_locs(i_loc), wrf_lon, wrf_lat);
                                else
                                    xx = misc_emissions_analysis.find_indices_in_radius_around_loc(base_locs(i_loc), wrf_lon, wrf_lat);
                                end

                                if is_2d
                                    data_for_today{1, i_file, i_var, i_loc} = wrf_value(xx);
                                else
                                    if read_upwind
                                        this_wrf_value = wrf_value(:,xx);
                                        if  isempty(this_wrf_value)
                                            data_for_today{1, i_file, i_var, i_loc} = wrf_value(:,xx);
                                        else
                                            indx = nansum(this_wrf_value(1:5,:),1) == min(nansum(this_wrf_value(1:5,:),1));
                                            data_for_today{1, i_file, i_var, i_loc} = this_wrf_value(:, indx);
                                        end
                                    else
                                        data_for_today{1, i_file, i_var, i_loc} = wrf_value(:,xx);
                                    end

                                end
                            end
                        end

                    end
                end
                
                fprintf('Average to locations...\n');
                for i_time =1:n_times
                    xx_date = winds_data{i_time}.dvec == this_date;
                    if sum(xx_date) < 1
                        fprintf('  %s not in time period %d\n', datestr(this_date), i_time);
                        continue
                    elseif sum(xx_date) > 1
                        E.notimplemented('Date matched multiple dates in location file')
                    end
                    read_qhcho = true;
                    for i_loc=1:n_locs
                        xx_hours = winds_data{i_time}.locs(i_loc).WindUsedBool(xx_date, :);
                        if ~skip_wrf
                            if ~isempty(winds_data{1}.locs(1).WindDir(xx_date,xx_hours))
                                winds_data{i_time}.locs(i_loc).WRFData.(field_2d).WindDir{xx_date} = nanmean(winds_data{1}.locs(1).WindDir(xx_date,xx_hours));
                                winds_data{i_time}.locs(i_loc).WRFData.(field_2d).WindSpeed{xx_date} = nanmean(winds_data{1}.locs(1).WindSpeed(xx_date,xx_hours));
                            else
                                winds_data{i_time}.locs(i_loc).WRFData.(field_2d).WindDir{xx_date} = nan;
                                winds_data{i_time}.locs(i_loc).WRFData.(field_2d).WindSpeed{xx_date} = nan;
                            end
                        end
                        for i_var=1:n_vars
                            if exist_vars_flag(i_time, i_var)
                                %fprintf('The variable is already exist%s\n', wrf_vars{i_var});
                                continue;
                            end
                            if any(regexp(wrf_vars{i_var}, '_vcd$'))
                                this_species = strrep(wrf_vars{i_var}, '_vcd','');
                                % todo: can't handle multiple years
                                profile = winds_data{1}.locs(i_loc).WRFData.Profile.(this_species){xx_date}*1e-6;
                                pres = winds_data{1}.locs(i_loc).WRFData.Profile.pres{xx_date};
                                vcd = nan(1,size(profile, 2));
                                for i_grid = 1:size(profile, 2)
                                    if all(isnan(profile(:,i_grid)))
                                        vcd(i_grid)= nan;
                                    else
                                        vcd(i_grid) = integPr2(profile(:,i_grid), pres(:,i_grid), pres(1,i_grid));
                                    end
                                end
                                 winds_data{i_time}.locs(i_loc).WRFData.(field_2d).(wrf_vars{i_var}){ xx_date} = vcd;
                                continue;
                            end
                            
                            if any(regexp(wrf_vars{i_var}, '\w+[^_vcd]\_aks$'))
                                this_species = strrep(wrf_vars{i_var}, '_aks', '');
                                if strcmpi(this_species, 'hcho')
                                    
                                    xx = misc_emissions_analysis.find_indices_in_radius_around_loc(winds_data{i_time}.locs(i_loc), wrf_lon, wrf_lat);
                                    site_lons = wrf_lon(xx);
                                    site_lats = wrf_lat(xx);
                                    try
                                        [closest_aks, closest_pres] = misc_ml_oh_prepare.match_closest_hcho_ak(this_date, site_lons, site_lats);

                                        if isnan(closest_pres)
                                            winds_data{i_time}.locs(i_loc).WRFData.(field_3d).(wrf_vars{i_var}){ xx_date} = nan;
                                        else
                                            for i_grid = 1:numel(site_lons)
                                                this_closest_pres = closest_pres(:,i_grid);
                                                this_closest_aks = closest_aks(:,i_grid);

                                                xx = ~isnan(this_closest_pres) & ~isnan(this_closest_aks);
                                                if numel(this_closest_pres(xx)) >=2 && ~isempty(winds_data{i_time}.locs(i_loc).WRFData.(field_3d).pres{ xx_date})
                                                    wrf_pres =  winds_data{i_time}.locs(i_loc).WRFData.(field_3d).pres{ xx_date}(:,i_grid); 
                                                    this_intrp_aks = interp1(this_closest_pres(xx), this_closest_aks(xx), wrf_pres,'nearest','extrap');
                                                    this_intrp_aks(this_intrp_aks<=0) = nan;
                                                    winds_data{i_time}.locs(i_loc).WRFData.(field_3d).(wrf_vars{i_var}){ xx_date}(:,i_grid) = ...
                                                        this_intrp_aks;
                                                else
                                                    winds_data{i_time}.locs(i_loc).WRFData.(field_3d).(wrf_vars{i_var}){ xx_date}(:,i_grid) = nan;
                                                end
                                            end
                                        end
                                    catch err
                                        winds_data{i_time}.locs(i_loc).WRFData.(field_3d).(wrf_vars{i_var}){ xx_date} = nan;
                                    end
                                    continue;
                                elseif strcmpi(this_species, 'qhcho')
                                    if read_qhcho
                                        fprintf('Reading QA4ECV HCHO column\n');
                                        qhcho_data = misc_ml_oh_prepare.read_qa4ecv_hcho_pixels(this_date);
                                        read_qhcho = false;
                                    end
                                    xx = misc_emissions_analysis.find_indices_in_radius_around_loc(winds_data{i_time}.locs(i_loc), wrf_lon, wrf_lat);
                                    site_lons = wrf_lon(xx);
                                    site_lats = wrf_lat(xx);
                                    try
                                        [closest_aks, closest_pres] = misc_ml_oh_prepare.match_closest_qa4ecv_hcho_ak_from_data(qhcho_data, site_lons, site_lats);

                                        if isnan(closest_pres)
                                            winds_data{i_time}.locs(i_loc).WRFData.(field_3d).(wrf_vars{i_var}){ xx_date} = nan;
                                        else
                                            for i_grid = 1:numel(site_lons)
                                                this_closest_pres = closest_pres(:,i_grid);
                                                this_closest_aks = closest_aks(:,i_grid);

                                                xx = ~isnan(this_closest_pres) & ~isnan(this_closest_aks);
                                                if numel(this_closest_pres(xx)) >=2 && ~isempty(winds_data{i_time}.locs(i_loc).WRFData.(field_3d).pres{ xx_date})
                                                    wrf_pres =  winds_data{i_time}.locs(i_loc).WRFData.(field_3d).pres{ xx_date}(:,i_grid); 
                                                    this_intrp_aks = interp1(this_closest_pres(xx), this_closest_aks(xx), wrf_pres,'nearest','extrap');
                                                    this_intrp_aks(this_intrp_aks<=0) = nan;
                                                    winds_data{i_time}.locs(i_loc).WRFData.(field_3d).(wrf_vars{i_var}){ xx_date}(:,i_grid) = ...
                                                        this_intrp_aks;
                                                else
                                                    wrf_pres =  winds_data{i_time}.locs(i_loc).WRFData.(field_3d).pres{ xx_date}(:,i_grid); 
                                                    winds_data{i_time}.locs(i_loc).WRFData.(field_3d).(wrf_vars{i_var}){ xx_date}(:,i_grid) = nan(size(wrf_pres));
                                                end
                                            end
                                        end
                                    catch err
                                        winds_data{i_time}.locs(i_loc).WRFData.(field_3d).(wrf_vars{i_var}){ xx_date} = nan;
                                    end
                                    continue;
                                elseif strcmpi(this_species, 'no2')
                                    [closest_aks, closest_pres] = misc_ml_oh_prepare.find_closest_behr_no2_ak(this_date, winds_data{i_time}.locs(i_loc).Longitude, ...
                                                                                                        winds_data{i_time}.locs(i_loc).Latitude);
                                    wrf_pres =  winds_data{i_time}.locs(i_loc).WRFData.(field_3d).pres(:, xx_date);     
                                    xx = ~isnan(closest_pres) & ~isnan(closest_aks);
                                    if numel(closest_pres(xx)) >=2
                                        winds_data{i_time}.locs(i_loc).WRFData.(field_3d).(wrf_vars{i_var})(:, xx_date) = ...
                                            interp1(closest_pres(xx), closest_aks(xx), wrf_pres,'linear','extrap');
                                    else
                                        winds_data{i_time}.locs(i_loc).WRFData.(field_3d).(wrf_vars{i_var})(:, xx_date) = nan;
                                    end
                                    continue;
                                else
                                    E.notimplemented('Aks for this species is not implemented');
                                end
                            end
                            
                            if any(regexp(wrf_vars{i_var}, '_vcd_aks$'))
                                this_species = strrep(wrf_vars{i_var}, '_vcd_aks','');
                                if any(strfind(this_species, '_surr'))
                                    this_species_aks = strcat(strrep(this_species,'_surr',''), '_aks');
                                elseif any(strfind(this_species, '_upwind'))
                                    this_species_aks = strcat(strrep(this_species,'_upwind',''), '_aks');
                                else
                                    this_species_aks = strcat(this_species, '_aks');
                                end
                                
                                if strcmpi(this_species,'qhcho')
                                    this_species = 'hcho';
                                end
                                % todo: can't handle multiple years
                                if isempty(winds_data{1}.locs(i_loc).WRFData.Profile.(this_species){xx_date})
                                     winds_data{i_time}.locs(i_loc).WRFData.(field_2d).(wrf_vars{i_var}){ xx_date} = nan;
                                else
                                    profile = winds_data{1}.locs(i_loc).WRFData.Profile.(this_species){xx_date}*1e-6 ...
                                        .* winds_data{1}.locs(i_loc).WRFData.Profile.(this_species_aks){xx_date} ;
                                    pres = winds_data{1}.locs(i_loc).WRFData.Profile.pres{xx_date};
                                    if all(isnan(profile))
                                        winds_data{i_time}.locs(i_loc).WRFData.(field_2d).(wrf_vars{i_var}){ xx_date} = nan;
                                    else
                                        this_grid = nan(1, size(pres,2));
                                        for i_grid = 1:size(pres,2)
                                            if all(isnan(profile(1:5,i_grid)))
                                                this_grid(i_grid) = nan;
                                            else
                                                this_grid(i_grid) = integPr2(profile(:,i_grid), pres(:,i_grid), pres(1,i_grid));
                                            end
                                        end
                                        winds_data{i_time}.locs(i_loc).WRFData.(field_2d).(wrf_vars{i_var}){ xx_date} = this_grid;
                                    end
                                end
                                continue;
                            end
                            
                            is_2d = is_var_2d(i_var);
                            if is_2d
                                loc_var_value = nanmean(cat(3, data_for_today{1, xx_hours, i_var, i_loc}),3);
                                substruct = field_2d;
                            else
                                loc_var_value = nanmean(cat(3, data_for_today{1, xx_hours, i_var, i_loc}),3);
                                substruct = field_3d;
                            end
                            winds_data{i_time}.locs(i_loc).WRFData.(substruct).(wrf_vars{i_var}){xx_date} = loc_var_value;
                        end
                    end
                end
            end
            
            for i_time=1:n_times
                savename = winds_data{i_time}.savename;
                savedata = rmfield(winds_data{i_time}, 'savename');
                save(savename, '-v7.3', '-struct', 'savedata' );
            end

            function yn = is_var_2d(ind)
                yn = ind >= n_3d_vars + 1;
            end
        end
        
        function make_location_small_radius_wrf_avgs_increment_vcd_file(time_periods,  varargin)
            E = JLLErrors;
            p = advInputParser;
            p.addParameter('overwrite_vars', {});
            p.addParameter('loc_inds', 1:49);
            p.parse(varargin{:});
            pout = p.Results;
            overwrite_vars = pout.overwrite_vars;
            
            loc_inds = pout.loc_inds;
            vcd_vars = {'behr_no2',  'qa4ecv_hcho', 'nasa_hcho', 'qa4ecv_l3_hcho'};
            base_locs = misc_emissions_analysis.read_locs_file();
            base_locs = misc_emissions_analysis.cutdown_locs_by_index(base_locs, loc_inds);
            base_locs_names = {base_locs.ShortName};
            
            wrf_data = load(misc_wrf_lifetime_analysis.wrf_avg_name(2014,'no2'));
            wrf_xlon = double(wrf_data.xlon);
            wrf_xlat = double(wrf_data.xlat);
            
            n_vars = numel(vcd_vars);
            n_locs = numel(base_locs);
            n_times = numel(time_periods);
            exist_vars_flag = false(n_times, n_vars);
            winds_data = cell(n_times,1);
            for i_time = 1:n_times
                [start_date, end_date] = misc_emissions_analysis.select_start_end_dates(time_periods{i_time});
                if isfile(misc_emissions_analysis.omi_small_radius_data_file_name(start_date, end_date))
                    winds_data{i_time} = load(misc_emissions_analysis.omi_small_radius_data_file_name(start_date, end_date));
                    increment = true;
                else
                    if isfile(misc_emissions_analysis.wrf_data_file_name(start_date, end_date))
                        winds_data{i_time} = load(misc_emissions_analysis.wrf_data_file_name(start_date, end_date));
                    else
                        winds_data{i_time} = load(misc_emissions_analysis.winds_file_name(start_date, end_date));
                    end
                    increment = false;
                end
                if exist('loc_inds', 'var')
                    winds_data{i_time}.locs = misc_emissions_analysis.cutdown_locs_by_index(winds_data{i_time}.locs, loc_inds);
                end
                % Double check that the winds structs match the base locs
                tmp_locs = winds_data{i_time}.locs;
                loc_names = {tmp_locs.ShortName};
                if ~isequal(loc_names, base_locs_names)
                    E.notimplemented('Winds locs different from base locs')
                end
                
                
                % Append a substructure to put the WRF data into. Have one
                % value per day, we'll average together whichever hours are
                % used.
                
                for i_loc=1:n_locs
                    for i_var = 1:numel(vcd_vars)
                        data_struct_2d.(vcd_vars{i_var}) = cell(1, size(tmp_locs(i_loc).WindUsedBool,1));
                    end
                    tmp_locs(i_loc).OMIData.VCD = data_struct_2d;
                    if increment
                        for i_var = 1:n_vars
                            exist_vars = fieldnames(winds_data{i_time}.locs(i_loc).OMIData.VCD);
                            var_name = vcd_vars{i_var};
                            if any(strcmpi(exist_vars, var_name))
                                if any(strcmpi(overwrite_vars, var_name))
                                    fprintf('    Overwrite %s\n', vcd_vars{i_var});
                                else
                                    fprintf('Reading %s from existing file.\n', vcd_vars{i_var});
                                    tmp_locs(i_loc).OMIData.VCD.(var_name)= winds_data{i_time}.locs(i_loc).OMIData.VCD.(var_name);
                                    exist_vars_flag(i_time, i_var) = true;
                                end
                            end
                        end
                    end
                end
                winds_data{i_time}.locs = tmp_locs;
                winds_data{i_time}.savename = misc_emissions_analysis.omi_small_radius_data_file_name(start_date, end_date);
            end
            
            all_years = unique(veccat(time_periods{:}));
            [total_starts, total_ends] = misc_emissions_analysis.select_start_end_dates(all_years);
            total_dvec = misc_emissions_analysis.make_datevec(total_starts, total_ends);

            common_opts = {'DEBUG_LEVEL', 1, 'dayofweek', 'UMTWRFS'};
            for d=1:numel(total_dvec)
                this_date = total_dvec(d);
                fprintf('Now on %s (%d of %d)\n', datestr(this_date), d, numel(total_dvec));
                
                fprintf('Average to locations...\n');
                for i_time = 1:n_times
                    xx_date = winds_data{i_time}.dvec == this_date;
                    if sum(xx_date) < 1
                        fprintf('  %s not in time period %d\n', datestr(this_date), i_time);
                        continue
                    elseif sum(xx_date) > 1
                        E.notimplemented('Date matched multiple dates in location file')
                    end
                    for i_var=1:n_vars

                        if increment && exist_vars_flag(i_time, i_var)
                            fprintf('The variable is already exist%s\n', vcd_vars{i_var});
                            continue;
                        end
                        if any(strcmpi(vcd_vars{i_var}, {'behr_no2','behr_no2_surr'}))
                            try
                                [value, lon, lat, weights] = behr_time_average(this_date, this_date, 'prof_mode', 'daily', common_opts{:});
                                no2.value = reshape(griddata(lon(:), lat(:), value(:), wrf_xlon(:), wrf_xlat(:)), size(wrf_xlon));
                                skip_behr_no2 = false;
                            catch err
                                fprintf('Error in behr no2, skip this date %s.\n', datestr(this_date));
                                skip_behr_no2 = true;
                            end
                        elseif any(strcmpi(vcd_vars{i_var}, {'omi_hcho', 'omi_hcho_surr'}))
                            try
                                [hcho.value, hcho.lon, hcho.lat] = omhcho_time_average(this_date, this_date, common_opts{:});
                                skip_omi_hcho = false;
                            catch err
                                fprintf('Error in omi hcho, skip this date %s.\n', datestr(this_date));
                                skip_omi_hcho = true;
                            end
            
                        elseif strcmpi(vcd_vars{i_var}, 'omi_o3')
                            try
                                [o3.value, o3.lon, o3.lat] = omo3_time_average(this_date, this_date, common_opts{:});
                                skip_omi_o3 = false;
                            catch err
                                fprintf('Error in omi o3, skip this date %s.\n', datestr(this_date));
                                skip_omi_o3 = true;
                            end
                        elseif strcmpi(vcd_vars{i_var}, 'qa4ecv_hcho')
                            try
                                [qhcho.value] = misc_ml_oh_prepare.qa4ecv_hcho_time_average_wrf_grids(this_date, this_date);
                                skip_qhcho = false;
                            catch err
                                fprintf('Error in qa4ecv hcho, skip this date %s.\n', datestr(this_date));
                                skip_qhcho = true;
                            end
                        elseif strcmpi(vcd_vars{i_var}, 'nasa_hcho')
                            try
                                [nhcho.value] = misc_ml_oh_prepare.nasa_hcho_time_average_wrf_grids(this_date, this_date);
                                skip_nhcho = false;
                            catch err
                                fprintf('Error in qa4ecv hcho, skip this date %s.\n', datestr(this_date));
                                skip_nhcho = true;
                            end
                        elseif strcmpi(vcd_vars{i_var}, 'qa4ecv_l3_hcho')
                            try
                                [nhcho.value] = misc_ml_oh_prepare.qa4ecv_hcho_time_average_wrf_grids_v2(this_date, this_date);
                                skip_qhcho_l3 = false;
                            catch err
                                fprintf('Error in qa4ecv hcho, skip this date %s.\n', datestr(this_date));
                                skip_qhcho_l3 = true;
                            end        
                        end
                        for i_loc=1:n_locs
                            if strcmpi(vcd_vars{i_var}, 'behr_no2')
                                if skip_behr_no2
                                    winds_data{i_time}.locs(i_loc).OMIData.VCD.(vcd_vars{i_var}){ xx_date} = nan(size(wrf_xlon(xx_radius)));
                                else
                                    xx_radius = misc_emissions_analysis.find_indices_in_radius_around_loc(winds_data{i_time}.locs(i_loc), wrf_xlon, wrf_xlat);
                                    winds_data{i_time}.locs(i_loc).OMIData.VCD.(vcd_vars{i_var}){ xx_date} = no2.value(xx_radius);
                                end
                            elseif strcmpi(vcd_vars{i_var}, 'behr_no2_surr')
                                if skip_behr_no2
                                    winds_data{i_time}.locs(i_loc).OMIData.VCD.(vcd_vars{i_var})( xx_date) = nan;
                                else
                                    xx_radius = misc_emissions_analysis.find_indices_in_radius_surr_loc(winds_data{i_time}.locs(i_loc), no2.lon, wrf_xlon, wrf_xlat);
                                    winds_data{i_time}.locs(i_loc).OMIData.VCD.(vcd_vars{i_var})( xx_date) = nanmean(no2.value(xx_radius));
                                end
                            elseif strcmpi(vcd_vars{i_var}, 'omi_hcho')
                                if skip_omi_hcho || isempty(hcho.value)
                                    winds_data{i_time}.locs(i_loc).OMIData.VCD.(vcd_vars{i_var})( xx_date) = nan;
                                else
                                    xx_radius = misc_emissions_analysis.find_indices_in_radius_around_loc(winds_data{i_time}.locs(i_loc), wrf_xlon, wrf_xlat);
                                    winds_data{i_time}.locs(i_loc).OMIData.VCD.(vcd_vars{i_var})( xx_date) = nanmean(hcho.value(xx_radius));
                                end
                            elseif strcmpi(vcd_vars{i_var}, 'omi_hcho_surr')
                                if skip_omi_hcho || isempty(hcho.value)
                                    winds_data{i_time}.locs(i_loc).OMIData.VCD.(vcd_vars{i_var})( xx_date) = nan;
                                else
                                    xx_radius = misc_emissions_analysis.find_indices_in_radius_surr_loc(winds_data{i_time}.locs(i_loc), wrf_xlon, wrf_xlat);
                                    winds_data{i_time}.locs(i_loc).OMIData.VCD.(vcd_vars{i_var})( xx_date) = nanmean(hcho.value(xx_radius));
                                end
                            elseif strcmpi(vcd_vars{i_var}, 'omi_o3')
                                if skip_omi_o3 || isempty(o3.value)
                                    winds_data{i_time}.locs(i_loc).OMIData.VCD.(vcd_vars{i_var})( xx_date) = nan;
                                else
                                    xx_radius = misc_emissions_analysis.find_indices_in_radius_around_loc(winds_data{i_time}.locs(i_loc), wrf_xlon, wrf_xlat);
                                    winds_data{i_time}.locs(i_loc).OMIData.VCD.(vcd_vars{i_var}){ xx_date} = o3.value(xx_radius);
                                end    
                            elseif strcmpi(vcd_vars{i_var}, 'qa4ecv_hcho')
                                if skip_qhcho
                                    winds_data{i_time}.locs(i_loc).OMIData.VCD.(vcd_vars{i_var}){ xx_date} =  nan(size(wrf_xlon(xx_radius)));
                                else
                                    xx_radius = misc_emissions_analysis.find_indices_in_radius_around_loc(winds_data{i_time}.locs(i_loc), wrf_xlon, wrf_xlat);
                                    winds_data{i_time}.locs(i_loc).OMIData.VCD.(vcd_vars{i_var}){ xx_date} = qhcho.value(xx_radius);
                                end
                            elseif strcmpi(vcd_vars{i_var}, 'nasa_hcho')
                                if skip_nhcho
                                    winds_data{i_time}.locs(i_loc).OMIData.VCD.(vcd_vars{i_var}){ xx_date} =  nan(size(wrf_xlon(xx_radius)));
                                else
                                    xx_radius = misc_emissions_analysis.find_indices_in_radius_around_loc(winds_data{i_time}.locs(i_loc), wrf_xlon, wrf_xlat);
                                    winds_data{i_time}.locs(i_loc).OMIData.VCD.(vcd_vars{i_var}){ xx_date} = nhcho.value(xx_radius);
                                end
                            elseif strcmpi(vcd_vars{i_var}, 'qa4ecv_l3_hcho')
                                if skip_qhcho_l3
                                    winds_data{i_time}.locs(i_loc).OMIData.VCD.(vcd_vars{i_var}){ xx_date} =  nan(size(wrf_xlon(xx_radius)));
                                else
                                    xx_radius = misc_emissions_analysis.find_indices_in_radius_around_loc(winds_data{i_time}.locs(i_loc), wrf_xlon, wrf_xlat);
                                    winds_data{i_time}.locs(i_loc).OMIData.VCD.(vcd_vars{i_var}){ xx_date} = nhcho.value(xx_radius);
                                end
                            end
                        end
                    end
                end
            end
            
            for i_time=1:n_times
                savename = winds_data{i_time}.savename;
                savedata = rmfield(winds_data{i_time}, 'savename');
                save(savename, '-v7.3', '-struct', 'savedata' );
            end

            function yn = is_var_2d(ind)
                yn = ind <= n_2d_vars;
            end
        end
        
        function merge_smallradius_omi_files(time_periods)
            E = JLLErrors;
            if ~exist('DEBUG_LEVEL', 'var')
                DEBUG_LEVEL = 2;
            end
            n_times = numel(time_periods);
            dvec_all = [];
            create_array = true;
            for i_time = 1:n_times
                [start_date, end_date] = misc_emissions_analysis.select_start_end_dates(time_periods{i_time});
                wrf_data{i_time} = load(misc_emissions_analysis.omi_small_radius_data_file_name(start_date, end_date));
                dvec_all = cat(2, dvec_all, wrf_data{i_time}.dvec);
                fprintf('Merging %s\n',misc_emissions_analysis.omi_small_radius_data_file_name(start_date, end_date));
                if create_array
                    loc_all = wrf_data{i_time}.locs;
                    create_array = false;
                else
                    omi_vars = fieldnames(loc_all(1).OMIData.VCD);
                    for i_loc = 1:numel(loc_all)
                        
                        for i_var=1:numel(omi_vars)
                            prev = loc_all(i_loc).OMIData.VCD.(omi_vars{i_var});
                            this_time = wrf_data{i_time}.locs(i_loc).OMIData.VCD.(omi_vars{i_var});
                            try 
                                if isnan(this_time{:})
                                    this_time = nan(this_size);
                                end
                            catch err
                                this_size = size(this_time);
                            end
                            loc_all(i_loc).OMIData.VCD.(omi_vars{i_var}) = cat(2,prev, this_time);
                        end
                    end
                    
                end
            end
            
            % Combine the line density structures then clear out the
            % original to save memory (can be 10+ GB).
            LD_all.locs = loc_all;
            % We already checked that all the date vectors are the same
            LD_all.dvec = dvec_all;
            % Keep all the write dates
            LD_all.write_date = datestr(now);
            
            new_save_name = fullfile(misc_emissions_analysis.site_info_dir, 'site_omi_small_radius_data_merge');
            
            if exist(new_save_name, 'file')
                if ~ask_yn(sprintf('%s exists. Overwrite? ', new_save_name))
                    return
                end
            end
            fprintf('Saving merged file as %s\n', new_save_name);
            save(new_save_name, '-v7.3', '-struct', 'LD_all');
            
        end
        
        function merge_smallradius_wrf_files(time_periods)
            E = JLLErrors;
            if ~exist('DEBUG_LEVEL', 'var')
                DEBUG_LEVEL = 2;
            end
            n_times = numel(time_periods);
            dvec_all = [];
            create_array = true;
            for i_time = 1:n_times
                [start_date, end_date] = misc_emissions_analysis.select_start_end_dates(time_periods{i_time});
                wrf_data{i_time} = load(misc_emissions_analysis.wrf_small_radius_data_file_name(start_date, end_date));
                dvec_all = cat(2, dvec_all, wrf_data{i_time}.dvec);
                fprintf('Merging %s\n',misc_emissions_analysis.wrf_small_radius_data_file_name(start_date, end_date));
                if create_array
                    
                    loc_all = wrf_data{i_time}.locs;
                    create_array = false;
                    
                else
                    wrf_vars = fieldnames(loc_all(1).WRFData.Averaged);
                    for i_loc = 1:numel(loc_all)
                        
                        for i_var=1:numel(wrf_vars)
                            prev = loc_all(i_loc).WRFData.Averaged.(wrf_vars{i_var});
                            this_time = wrf_data{i_time}.locs(i_loc).WRFData.Averaged.(wrf_vars{i_var});
                            loc_all(i_loc).WRFData.Averaged.(wrf_vars{i_var}) = cat(2,prev, this_time);
                        end
                    end
                    
                    wrf_vars = fieldnames(loc_all(1).WRFData.Profile);
                    for i_loc = 1:numel(loc_all)
                        
                        for i_var=1:numel(wrf_vars)
                            prev = loc_all(i_loc).WRFData.Profile.(wrf_vars{i_var});
                            this_time = wrf_data{i_time}.locs(i_loc).WRFData.Profile.(wrf_vars{i_var});
                            loc_all(i_loc).WRFData.Profile.(wrf_vars{i_var}) = cat(2,prev, this_time);
                        end
                    end
                end
            end
            
            % Combine the line density structures then clear out the
            % original to save memory (can be 10+ GB).
            LD_all.locs = loc_all;
            % We already checked that all the date vectors are the same
            LD_all.dvec = dvec_all;
            % Keep all the write dates
            LD_all.write_date = datestr(now);
            
            new_save_name = fullfile(misc_emissions_analysis.site_info_dir, 'site_wrf_small_radius_data_merge');
            
            if exist(new_save_name, 'file')
                if ~ask_yn(sprintf('%s exists. Overwrite? ', new_save_name))
                    return
                end
            end
            fprintf('Saving merged file as %s\n', new_save_name);
            save(new_save_name, '-v7.3', '-struct', 'LD_all');
            
        end
        
        function prepare_ml_small_radius_output_nomerge(time_periods, loc_inds)
            E = JLLErrors;
            if ~exist('DEBUG_LEVEL', 'var')
                DEBUG_LEVEL = 2;
            end
            n_times = numel(time_periods);
            dvec_all = [];
            create_array = true;
            for i_time = 1:n_times
                [start_date, end_date] = misc_emissions_analysis.select_start_end_dates(time_periods{i_time});
                wrf_data{i_time} = load(misc_emissions_analysis.wrf_small_radius_data_file_name(start_date, end_date));
                dvec_all = cat(2, dvec_all, wrf_data{i_time}.dvec);
                fprintf('Merging %s\n',misc_emissions_analysis.wrf_small_radius_data_file_name(start_date, end_date));
                if create_array
                    
                    loc_all = wrf_data{i_time}.locs;
                    create_array = false;
                    
                else
                    wrf_vars = fieldnames(loc_all(1).WRFData.Averaged);
                    for i_loc = 1:numel(loc_all)
                        
                        for i_var=1:numel(wrf_vars)
                            prev = loc_all(i_loc).WRFData.Averaged.(wrf_vars{i_var});
                            this_time = wrf_data{i_time}.locs(i_loc).WRFData.Averaged.(wrf_vars{i_var});
                            loc_all(i_loc).WRFData.Averaged.(wrf_vars{i_var}) = cat(2,prev, this_time);
                        end
                    end
                    
                    wrf_vars = fieldnames(loc_all(1).WRFData.Profile);
                    for i_loc = 1:numel(loc_all)
                        
                        for i_var=1:numel(wrf_vars)
                            prev = loc_all(i_loc).WRFData.Profile.(wrf_vars{i_var});
                            this_time = wrf_data{i_time}.locs(i_loc).WRFData.Profile.(wrf_vars{i_var});
                            loc_all(i_loc).WRFData.Profile.(wrf_vars{i_var}) = cat(2,prev, this_time);
                        end
                    end
                end
            end
            

            base_locs = misc_ml_oh_prepare.make_wrf_lon_lat_grids(loc_inds);
            
            for i_loc = 1:numel(loc_inds)
                loc_ind = loc_inds(i_loc);
                loc = loc_all(loc_ind);
                fds = fieldnames(loc.WRFData.Averaged);
                create_dvec = true;
                geos_lon = base_locs(loc_ind).wrf_lon;
                geos_lat = base_locs(loc_ind).wrf_lat;
                for i=1:numel(fds)
                    this_data = [];
                    if create_dvec
                        this_dvec = [];
                        this_geos_lon = [];
                        this_geos_lat = [];
                    end
                    for i_day = 1:numel(loc.WRFData.Averaged.(fds{i}))
                        this_res = loc.WRFData.Averaged.(fds{i})(i_day);
                        snum(i_day) = numel(this_res{:});
                        if create_dvec
                            this_data = cat(1, this_data, this_res{:}(:));
                        else
                            if snum(i_day) ~= snum_save(i_day)
                                this_data = cat(1, this_data, nan(snum_save(i_day),1));
                            else
                                this_data = cat(1, this_data, this_res{:}(:));
                            end
                        end
                        
                        
                        if create_dvec
                            this_dvec = cat(1, this_dvec, zeros(size(this_res{:}(:)))+dvec_all(i_day));
                            if isempty(this_res{:}(:)) || all(isnan(this_res{:}(:)))
                                this_geos_lon = cat(1, this_geos_lon, []);
                                this_geos_lat = cat(1, this_geos_lat, []);
                            else
                                this_geos_lon = cat(1, this_geos_lon, geos_lon);
                                this_geos_lat = cat(1, this_geos_lat, geos_lat);
                            end
                        end
                    end
                    loc.WRFData.Averaged.(fds{i}) = this_data;
                    if create_dvec
                        loc.WRFData.Averaged.dvec = this_dvec;
                        loc.WRFData.Averaged.xlon = this_geos_lon;
                        loc.WRFData.Averaged.xlat = this_geos_lat;
                        snum_save = snum;
                        create_dvec = false;
                    end
                end
                loc.WRFData.Averaged = rmfield(loc.WRFData.Averaged, {'WindDir', 'WindSpeed', 'hcho_vcd_aks'});
                tab = struct2table(loc.WRFData.Averaged,'AsArray',false);
                save_name = sprintf('wrf_data_merged_%02d.csv', loc_ind);
                save_name = fullfile(misc_emissions_analysis.site_info_dir, 'ML-WRF-Small-Radius', save_name);
                writetable(tab,save_name,'Delimiter',',');
            end
            
        end
        
        function prepare_omi_ml_small_radius_output(loc_inds)
            data = load(fullfile(misc_emissions_analysis.site_info_dir, 'site_omi_small_radius_data_merge'));
            dvec = data.dvec;
            base_locs = misc_ml_oh_prepare.make_wrf_lon_lat_grids(loc_inds);
            
            for i_loc = 1:numel(loc_inds)
                loc_ind = loc_inds(i_loc);
                loc = data.locs(loc_ind);
                fds = fieldnames(loc.OMIData.VCD);
                create_dvec = true;
                geos_lon = base_locs(loc_ind).wrf_lon;
                geos_lat = base_locs(loc_ind).wrf_lat;
                size_read = loc.OMIData.VCD.behr_no2(1);
                this_day_size = numel(geos_lon);
                for i=1:numel(fds)
                    this_data = [];
                    if create_dvec
                        this_dvec = [];
                        this_geos_lon = [];
                        this_geos_lat = [];
                    end
                    
                    for i_day = 1:numel(loc.OMIData.VCD.(fds{i}))
                        this_res = loc.OMIData.VCD.(fds{i})(i_day);
                        if numel(this_res{:}) ~= this_day_size
                            this_res = {nan(this_day_size,1)};
                        end
                        this_data = cat(1, this_data, this_res{:}(:));
                        this_num(i, i_day) = numel(this_res{:}(:));
                        if create_dvec
                            this_dvec = cat(1, this_dvec, zeros(size(this_res{:}(:)))+dvec(i_day));
                            if isempty(this_res{:}(:))
                                this_geos_lon = cat(1, this_geos_lon, []);
                                this_geos_lat = cat(1, this_geos_lat, []);
                            else
                                this_geos_lon = cat(1, this_geos_lon, geos_lon);
                                this_geos_lat = cat(1, this_geos_lat, geos_lat);
                            end
                        end
                    end
                    loc.OMIData.VCD.(fds{i}) = this_data;
                    if create_dvec
                        loc.OMIData.VCD.dvec = this_dvec;
                        loc.OMIData.VCD.xlon = this_geos_lon;
                        loc.OMIData.VCD.xlat = this_geos_lat;
                        create_dvec = false;
                    end
                end
                tab = struct2table(loc.OMIData.VCD,'AsArray',false);
                save_name = sprintf('omi_data_merged_%02d.csv', loc_ind);
                save_name = fullfile(misc_emissions_analysis.site_info_dir, 'ML-OMI-Small-Radius', save_name);
                writetable(tab,save_name,'Delimiter',',');
            end
            
        end
        
        %%%%%%UTILITY
        function [closest_aks, closest_pres] = find_closest_behr_no2_ak(start_date, wrf_lon, wrf_lat)
            behr_dir = behr_paths.BEHRMatSubdir('us', 'daily');
            file_pattern = behr_filename(start_date, 'daily', 'us', 'mat');
            F = dir(fullfile(behr_dir, file_pattern));
            D = load(fullfile(behr_dir, F(1).name), 'Data');
            for i_swath = 1:numel(D.Data)
                lon = D.Data(i_swath).Longitude;
                lat = D.Data(i_swath).Latitude;
                is_in = lon >= wrf_lon - 1 & lon <= wrf_lon + 1 & lat >= wrf_lat - 1 & lat <= wrf_lat + 1;
                xx = any(is_in, 1);
                if ~any(xx)
                    continue
                end
                aks = D.Data(i_swath).BEHRAvgKernels;
                pres = D.Data(i_swath).BEHRPressureLevels;
                indx = (lon - wrf_lon).^2 + (lat-wrf_lat).^2 ==min((lon(:) - wrf_lon).^2 + (lat(:)-wrf_lat).^2);
                closest_aks = nanmean(aks( :, indx),2);
                closest_pres = nanmean(pres( :,indx),2);
                return
            end
            closest_aks = nan;
            closest_pres = nan;
        end
        
        function [closest_aks, closest_pres] = find_closest_hcho_ak(start_date, wrf_lon, wrf_lat)
            dnum =  validate_date(start_date);
            last_path = '';
            curr_files = get_file_for_date(dnum);
            for i_file = 1:numel(curr_files)
                this_file = fullfile(last_path, curr_files(i_file).name);
                xx_rows = in_domain(this_file);
                if ~any(xx_rows)
                    continue
                end
                Data = load_omhcho_h5(this_file, xx_rows);
                for i_vel = 1:size(Data.ScatteringWeights, 3)
                    aks(:,:,i_vel) = Data.ScatteringWeights(:,:,i_vel)./Data.AirMassFactor;
                end
                aks = permute(aks,[3 1 2]);
                pres = Data.ClimatologyLevels;
                pres = permute(pres, [3 1 2]);
                pixel_lon = Data.Longitude;
                pixel_lat = Data.Latitude;
                indx = (pixel_lon - wrf_lon).^2 + (pixel_lat-wrf_lat).^2 ==min((pixel_lon(:) - wrf_lon).^2 + (pixel_lat(:)-wrf_lat).^2);
                closest_aks = nanmean(aks( :, indx),2);
                closest_pres = nanmean(pres( :,indx),2);
                return
            end
            closest_aks = nan;
            closest_pres = nan;
            
            function F = get_file_for_date(dnum)
                year_str = datestr(dnum, 'yyyy');
                month_str = datestr(dnum, 'mm');
                curr_path = fullfile(behr_paths.omno2_dir, '..', '..', 'OMHCHO', year_str, month_str);
                if ~strcmp(curr_path, last_path)
                    % Use cached list of files if we can to speed things up
                    omhcho_files = dir(fullfile(curr_path, '*.he5'));
                end
                last_path = curr_path;

                xx = regcmp({omhcho_files.name}, sprintf('^OMI-Aura_L2-OMHCHO_%sm%s%s', year_str, month_str, datestr(dnum, 'dd')));
                F = omhcho_files(xx);
            end
    
            function xx = in_domain(h5filename)
                hi = h5info(h5filename);
                lon = h5read(hi.Filename, h5dsetname(hi,1,2,1,2,'Longitude'));
                lat = h5read(hi.Filename, h5dsetname(hi,1,2,1,2,'Latitude'));
                is_in = lon >= wrf_lon - 1 & lon <= wrf_lon + 1 & lat >= wrf_lat - 1 & lat <= wrf_lat + 1;
                xx = any(is_in, 1);
            end

            function Data = load_omhcho_h5(h5filename, xx_rows)
                hi = h5info(h5filename);
                group1_fields = {'AirMassFactor'};
                group2_fields = {'Longitude', 'Latitude'};
                group3_fields = {'ScatteringWeights', 'ClimatologyLevels'};
                Data = make_empty_struct_from_cell([group1_fields, group2_fields, group3_fields]);

                xx_corner_rows = false(1,length(xx_rows)+1);
                xx_corner_rows(xx_rows) = true;
                xx_corner_rows(find(xx_rows,1,'last')+1) = true;

                for i = 1:numel(group1_fields)
                    fname = group1_fields{i};
                    Data.(fname) = h5readomi(hi.Filename, h5dsetname(hi,1,2,1,1, fname));
                    if size(Data.(fname),2) == length(xx_rows)
                        Data.(fname)(:,~xx_rows) = [];
                    elseif size(Data.(fname),2) == length(xx_corner_rows)
                        Data.(fname)(:,~xx_corner_rows) = [];
                    else
                        E.notimplemented('size(dataset,2) = %d', size(Data.(fname),2));
                    end
                end

                for i = 1:numel(group2_fields)
                    fname = group2_fields{i};
                    Data.(fname) = h5readomi(hi.Filename, h5dsetname(hi,1,2,1,2, fname), 'keep_type', true);
                    Data.(fname)(:,~xx_rows) = [];
                end

                for i = 1:numel(group3_fields)
                    fname = group3_fields{i};
                    Data.(fname) = h5readomi(hi.Filename, h5dsetname(hi,1,2,1,1, fname));
                    if size(Data.(fname),2) == length(xx_rows)
                        Data.(fname)(:,~xx_rows,:) = [];
                    elseif size(Data.(fname),2) == length(xx_corner_rows)
                        Data.(fname)(:,~xx_corner_rows,:) = [];
                    else
                        E.notimplemented('size(dataset,2) = %d', size(Data.(fname),2));
                    end
                end
            end
        end
        
        function [closest_aks, closest_pres] = match_closest_hcho_ak(start_date, site_lons, site_lats)
            dnum =  validate_date(start_date);
            last_path = '';
            curr_files = get_file_for_date(dnum);
            for i_file = 1:numel(curr_files)
                this_file = fullfile(last_path, curr_files(i_file).name);
                xx_rows = in_domain(this_file);
                if ~any(xx_rows)
                    continue
                end
                Data = load_omhcho_h5(this_file, xx_rows);
                for i_vel = 1:size(Data.ScatteringWeights, 3)
                    aks(:,:,i_vel) = Data.ScatteringWeights(:,:,i_vel)./Data.AirMassFactor;
                end
                aks = permute(aks,[3 1 2]);
                pres = Data.ClimatologyLevels;
                pres = permute(pres, [3 1 2]);
                pixel_lon = Data.Longitude;
                pixel_lat = Data.Latitude;
                
                for i_site = 1:numel(site_lons)
                    indx = (pixel_lon - site_lons(i_site)).^2 + (pixel_lat-site_lats(i_site)).^2 ==min((pixel_lon(:) - site_lons(i_site)).^2 + (pixel_lat(:)-site_lats(i_site)).^2);
                    closest_aks(:,i_site) = nanmean(aks( :, indx),2);
                    closest_pres(:,i_site)  = nanmean(pres( :,indx),2);
                end
                
                return
            end
            closest_aks = nan;
            closest_pres = nan;
            
            function F = get_file_for_date(dnum)
                year_str = datestr(dnum, 'yyyy');
                month_str = datestr(dnum, 'mm');
                curr_path = fullfile(behr_paths.omno2_dir, '..', '..', 'OMHCHO', year_str, month_str);
                if ~strcmp(curr_path, last_path)
                    % Use cached list of files if we can to speed things up
                    omhcho_files = dir(fullfile(curr_path, '*.he5'));
                end
                last_path = curr_path;

                xx = regcmp({omhcho_files.name}, sprintf('^OMI-Aura_L2-OMHCHO_%sm%s%s', year_str, month_str, datestr(dnum, 'dd')));
                F = omhcho_files(xx);
            end
    
            function xx = in_domain(h5filename)
                hi = h5info(h5filename);
                lon = h5read(hi.Filename, h5dsetname(hi,1,2,1,2,'Longitude'));
                lat = h5read(hi.Filename, h5dsetname(hi,1,2,1,2,'Latitude'));
                is_in = lon >= min(site_lons) - 1 & lon <= max(site_lons) + 1 & lat >= min(site_lats) - 1 & lat <= max(site_lats) + 1;
                xx = any(is_in, 1);
            end

            function Data = load_omhcho_h5(h5filename, xx_rows)
                hi = h5info(h5filename);
                group1_fields = {'AirMassFactor'};
                group2_fields = {'Longitude', 'Latitude'};
                group3_fields = {'ScatteringWeights', 'ClimatologyLevels'};
                Data = make_empty_struct_from_cell([group1_fields, group2_fields, group3_fields]);

                xx_corner_rows = false(1,length(xx_rows)+1);
                xx_corner_rows(xx_rows) = true;
                xx_corner_rows(find(xx_rows,1,'last')+1) = true;

                for i = 1:numel(group1_fields)
                    fname = group1_fields{i};
                    Data.(fname) = h5readomi(hi.Filename, h5dsetname(hi,1,2,1,1, fname));
                    if size(Data.(fname),2) == length(xx_rows)
                        Data.(fname)(:,~xx_rows) = [];
                    elseif size(Data.(fname),2) == length(xx_corner_rows)
                        Data.(fname)(:,~xx_corner_rows) = [];
                    else
                        E.notimplemented('size(dataset,2) = %d', size(Data.(fname),2));
                    end
                end

                for i = 1:numel(group2_fields)
                    fname = group2_fields{i};
                    Data.(fname) = h5readomi(hi.Filename, h5dsetname(hi,1,2,1,2, fname), 'keep_type', true);
                    Data.(fname)(:,~xx_rows) = [];
                end

                for i = 1:numel(group3_fields)
                    fname = group3_fields{i};
                    Data.(fname) = h5readomi(hi.Filename, h5dsetname(hi,1,2,1,1, fname));
                    if size(Data.(fname),2) == length(xx_rows)
                        Data.(fname)(:,~xx_rows,:) = [];
                    elseif size(Data.(fname),2) == length(xx_corner_rows)
                        Data.(fname)(:,~xx_corner_rows,:) = [];
                    else
                        E.notimplemented('size(dataset,2) = %d', size(Data.(fname),2));
                    end
                end
            end
        end
        
        function [closest_aks, closest_pres] = match_closest_qa4ecv_hcho_ak(start_date, site_lons, site_lats)
            dnum =  validate_date(start_date);
            last_path = '';
            curr_files = get_file_for_date(dnum);
            for i_file = 1:numel(curr_files)
                this_file = fullfile(last_path, curr_files(i_file).name);
                xx_rows = in_domain(this_file);
                if ~any(xx_rows)
                    continue
                end

                pixel_lon = ncread(this_file, '/PRODUCT/longitude');
                pixel_lat = ncread(this_file, '/PRODUCT/latitude');
                pres_a = ncread(this_file, '/PRODUCT/tm5_pressure_level_b');
                pres_b = ncread(this_file, '/PRODUCT/tm5_pressure_level_a');
                surf_pres = ncread(this_file, '/PRODUCT/tm5_surface_pressure');
                averaging_kernel = ncread(this_file, '/PRODUCT/averaging_kernel');
                
                for i_site = 1:numel(site_lons)
                    indx = (pixel_lon - site_lons(i_site)).^2 + (pixel_lat-site_lats(i_site)).^2 ==min((pixel_lon(:) - site_lons(i_site)).^2 + (pixel_lat(:)-site_lats(i_site)).^2);
                    closest_aks(:,i_site) = nanmean(averaging_kernel( :, indx),2);
                    closest_pres(:,i_site)  = nanmean(pres_a*nanmean(surf_pres(indx)) + pres_b/100, 1);
                end

                return
            end
            closest_aks = nan;
            closest_pres = nan;
            
            function F = get_file_for_date(dnum)
                year_str = datestr(dnum, 'yyyy');
                month_str = datestr(dnum, 'mm');
                curr_path = fullfile('/Volumes/share-sat/SAT/OMI/OMHCHO_QA4ECV_l2', year_str);
                if ~strcmp(curr_path, last_path)
                    % Use cached list of files if we can to speed things up
                    omhcho_files = dir(fullfile(curr_path, '*.nc'));
                end
                last_path = curr_path;
                xx = regcmp({omhcho_files.name}, sprintf('^QA4ECV_L2_HCHO_OMI_%s%s%s', year_str, month_str, datestr(dnum, 'dd')));
                F = omhcho_files(xx);
            end
    
            function xx = in_domain(ncfilename)
                hi = ncinfo(ncfilename);
                lon = ncread(hi.Filename, '/PRODUCT/longitude');
                lat = ncread(hi.Filename, '/PRODUCT/latitude');
                is_in = lon >= min(site_lons) - 1 & lon <= max(site_lons) + 1 & lat >= min(site_lats) - 1 & lat <= max(site_lats) + 1;
                xx = any(is_in, 1);
            end

        end
        
        function Data = read_qa4ecv_hcho_pixels(start_date)
            dnum =  validate_date(start_date);
            last_path = '';
            curr_files = get_file_for_date(dnum);
            pixel_lon = [];
            pixel_lat = [];
            pres = [];
            averaging_kernel = [];
            for i_file = 1:numel(curr_files)
                this_file = fullfile(last_path, curr_files(i_file).name);
                xx_rows = in_domain(this_file);
                if ~any(xx_rows)
                    continue
                end
                
                pixel_lon = cat(2, pixel_lon, ncread(this_file, '/PRODUCT/longitude'));
                pixel_lat = cat(2, pixel_lat, ncread(this_file, '/PRODUCT/latitude'));
                this_pres_a = nanmean(ncread(this_file, '/PRODUCT/tm5_pressure_level_a'),1)';
                this_pres_b = nanmean(ncread(this_file, '/PRODUCT/tm5_pressure_level_b'),1)';
                this_surf_pres = ncread(this_file, '/PRODUCT/tm5_surface_pressure');
                this_pres_a = repmat(this_pres_a, 1, size(this_surf_pres,1), size(this_surf_pres,2));
                this_pres_b = repmat(this_pres_b, 1, size(this_surf_pres,1), size(this_surf_pres,2));
                this_surf_pres = permute(repmat(this_surf_pres, 1,1,size(this_pres_a,1)),[3 1 2]);
                this_pres = this_pres_b.*this_surf_pres + this_pres_a/100;
                pres = cat(3, pres, this_pres);
                averaging_kernel =cat(3, averaging_kernel,  ncread(this_file, '/PRODUCT/averaging_kernel'));

            end
            Data.pixel_lon = pixel_lon;
            Data.pixel_lat = pixel_lat;
            Data.pres = pres;
            Data.averaging_kernel = averaging_kernel;
            
            function F = get_file_for_date(dnum)
                year_str = datestr(dnum, 'yyyy');
                month_str = datestr(dnum, 'mm');
                curr_path = fullfile('/Volumes/share-sat/SAT/OMI/OMHCHO_QA4ECV_l2', year_str);
                if ~strcmp(curr_path, last_path)
                    % Use cached list of files if we can to speed things up
                    omhcho_files = dir(fullfile(curr_path, '*.nc'));
                end
                last_path = curr_path;
                xx = regcmp({omhcho_files.name}, sprintf('^QA4ECV_L2_HCHO_OMI_%s%s%s', year_str, month_str, datestr(dnum, 'dd')));
                F = omhcho_files(xx);
            end
    
            function xx = in_domain(ncfilename)
                try
                    hi = ncinfo(ncfilename);
                    lon = ncread(hi.Filename, '/PRODUCT/longitude');
                    lat = ncread(hi.Filename, '/PRODUCT/latitude');
                    is_in = lon >= -125 & lon <= -65 & lat >= 25 & lat <= 50;
                    xx = any(is_in, 1);
                catch 
                    xx = 0;
                end
            end

        end
        
        function [closest_aks, closest_pres] = match_closest_qa4ecv_hcho_ak_from_data(Data, site_lons, site_lats)
            pixel_lon = Data.pixel_lon;
            pixel_lat = Data.pixel_lat;
            pres = Data.pres;
            averaging_kernel = Data.averaging_kernel;
            
            for i_site = 1:numel(site_lons)
                indx = (pixel_lon - site_lons(i_site)).^2 + (pixel_lat-site_lats(i_site)).^2 ==min((pixel_lon(:) - site_lons(i_site)).^2 + (pixel_lat(:)-site_lats(i_site)).^2);
                closest_aks(:,i_site) = nanmean(averaging_kernel( :, indx),2);
                closest_pres(:,i_site)  = nanmean(pres( :, indx),2);
            end

        end
        
        
    end
end