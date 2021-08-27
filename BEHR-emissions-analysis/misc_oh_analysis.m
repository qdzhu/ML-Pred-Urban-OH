classdef misc_oh_analysis
    %UNTITLED Summary of this class goes here
    %   Detailed explanation goes here
    
    properties(Constant)
        
    end
    
    methods(Static)
        %%%%%%%%%%%%%%%%%%%%%%%%%
        % Property-like methods %
        %%%%%%%%%%%%%%%%%%%%%%%%%
        
        function value = epa_data_dir
            value = fullfile(misc_emissions_analysis.workspace_dir, 'EPA-Data');
        end
        
        function windows = all_time_periods()
            wc = 2006:2013;
            windows = cell(size(wc));
            for i=1:numel(windows)
                windows{i} = (wc(i)-1):(wc(i)+1);
            end
        end
        
        function labels = all_time_period_labels()
            labels = cellfun(@sprintf_ranges, misc_oh_analysis.all_time_periods, 'uniform', false);
        end
        
        function db = rate_db()
            db_file = fullfile(misc_oh_analysis.epa_data_dir, 'voc_rates.sqlite');
            db = sqlite(db_file);
        end
        
        function [annual_files, years] = list_epa_files()
            annual_files = dirff(fullfile(misc_oh_analysis.epa_data_dir, 'annual*.csv'));
            years = regexp({annual_files.name}, '\d{4}(?=\.csv)', 'match', 'once');
            years = cellfun(@str2double, years);
        end
    
        %%%%%%%%%%%%%%%%%%%%
        % Plotting methods %
        %%%%%%%%%%%%%%%%%%%%
        
        function plot_oh_model_curves(varargin)
            p = advInputParser;
            p.addParameter('locs', misc_emissions_analysis.nine_cities);
            p.addParameter('oh_type', 'invert_hcho');
            
            p.parse(varargin{:});
            pout = p.Results;
            
            oh_type = pout.oh_type;
            loc_inds = misc_emissions_analysis.convert_input_loc_inds(pout.locs);
            
            n_locs = numel(loc_inds);
            n_times = numel(misc_oh_analysis.all_time_periods);
            locs = misc_emissions_analysis.read_locs_file();
            locs = misc_emissions_analysis.cutdown_locs_by_index(locs, loc_inds);
            
            % prepare 1 figure per location
            figs = gobjects(n_locs, 1);
            
            all_ax = gobjects(n_locs, 4);
            series = gobjects(n_locs, n_times, 2);
            ax_options = {'xscale', 'log'};
            for i_loc = 1:n_locs
                figs(i_loc) = figure;
                all_ax(i_loc,1) = subplot(2,2,1);
                set(all_ax(i_loc,1), ax_options{:});
                title(sprintf('%s weekdays', locs(i_loc).ShortName));
                xlabel('[NO_x] (ppb)'); ylabel('[OH] (ppt)');
                
                all_ax(i_loc,2) = subplot(2,2,2);
                set(all_ax(i_loc,2), ax_options{:});
                title(sprintf('%s weekends', locs(i_loc).ShortName));
                xlabel('[NO_x] (ppb)');
                
                % for P(HOx), VOCr, and alpha
                all_ax(i_loc,3) = subplot(2,2,3);
                all_ax(i_loc,4) = subplot(2,2,4);
                
                % resize to make top plots bigger
                all_ax(i_loc,1).Position(2) = 0.5;
                all_ax(i_loc,1).Position(4) = 0.4;
                all_ax(i_loc,2).Position(2) = 0.5;
                all_ax(i_loc,2).Position(4) = 0.4;
                all_ax(i_loc,3).Position(4) = 0.25;
                all_ax(i_loc,4).Position(4) = 0.25;
                
                subplot_stretch(1.5,2,'figh',figs(i_loc));
            end
            
            % Load each time period's OH file; get the VOCR, PHOx, alpha,
            % [OH] and [NOx]. Use the first three to make the OH vs NOx
            % steady state curve for that year, use [OH] and [NOx] to mark
            % on the plot where the analysis says that city is.
            year_labels = misc_oh_analysis.all_time_period_labels;
            phox_all = nan(n_locs, n_times, 2);
            vocr_all = nan(n_locs, n_times, 2);
            alpha_all = nan(n_locs, n_times, 2);
            
            for i_yr = 1:n_times
                fprintf('Working on %s\n', year_labels{i_yr});
                years = misc_oh_analysis.all_time_periods{i_yr};
                OH = load(misc_emissions_analysis.oh_file_name(years));
                year_color = misc_oh_analysis.get_time_period_color(years);
                for i_loc = 1:n_locs
                    ii_loc = loc_inds(i_loc);
                    wkday_ax = all_ax(i_loc, 1);
                    wkday_oh = OH.locs_wkday(ii_loc).OH.(oh_type);
                    [series(i_loc, i_yr, 1), phox_all(i_loc, i_yr, 1), alpha_all(i_loc, i_yr, 1), vocr_all(i_loc, i_yr, 1)]...
                        = plot_oh_curve(wkday_ax, wkday_oh, year_color);
                    
                    wkend_ax = all_ax(i_loc, 2);
                    wkend_oh = OH.locs_wkend(ii_loc).OH.(oh_type);
                    [series(i_loc, i_yr, 2), phox_all(i_loc, i_yr, 2), alpha_all(i_loc, i_yr, 2), vocr_all(i_loc, i_yr, 2)]...
                        = plot_oh_curve(wkend_ax, wkend_oh, year_color);
                end
            end
            
            % Go back and add the legend to the top two plots and plot the
            % model parameters (PHOx, alpha, VOCr) on the bottom plots
            for i_loc = 1:n_locs
                legend(all_ax(i_loc,1), series(i_loc, :, 1)', year_labels);
                legend(all_ax(i_loc,2), series(i_loc, :, 2)', year_labels);
                
                plot_model_params(all_ax(i_loc, 3), squeeze(phox_all(i_loc, :, 1)), squeeze(alpha_all(i_loc, :, 1)), squeeze(vocr_all(i_loc, :, 1)));
                plot_model_params(all_ax(i_loc, 4), squeeze(phox_all(i_loc, :, 2)), squeeze(alpha_all(i_loc, :, 2)), squeeze(vocr_all(i_loc, :, 2)));
            end
            
            
            
            function [l, phox, alpha, vocr] = plot_oh_curve(ax, loc_oh, color)
                line_opts = {'color', color, 'linewidth', 2};
                bad_vals = false;
                try
                    nox = loc_oh.nox;
                    oh = loc_oh.oh;
                    
                    phox = loc_oh.phox;
                    alpha = loc_oh.alpha;
                    vocr = loc_oh.vocr;
                catch err
                    if strcmp(err.identifier, 'MATLAB:nonExistentField')
                        fprintf('Have to skip this location/year: %s\n', err.message);
                        bad_vals = true;
                    else
                        rethrow(err)
                    end
                end
                
                if bad_vals || isnan(phox) || isnan(alpha) || isnan(vocr)
                    % create a dummy line to show up in the legend
                    l=line(ax, nan, nan, line_opts{:});
                    phox = nan;
                    alpha = nan;
                    vocr = nan;
                    return
                end
                
                x_nox = logspace(9,12,20);
                y_oh = nan(size(x_nox));
                fprintf('  solving HOx steady state: ');
                parfor i = 1:numel(x_nox)
                    fprintf('*');
                    y_oh(i) = hox_ss_solver(x_nox(i), phox, vocr, alpha);
                end
                fprintf('\n');
                
                l=line(ax, x_nox/2e10, y_oh/2e7, line_opts{:});
                line(ax, nox/2e10, oh/2e7, line_opts{:}, 'marker', 'o', 'markersize', 8, 'linestyle', 'none');
            end
            
            function plot_model_params(ax, loc_phox, loc_alpha, loc_vocr)
                x_years = cellfun(@mean, misc_oh_analysis.all_time_periods);
                x_labels = misc_oh_analysis.all_time_period_labels;
                
                l = gobjects(3,1);
                l(1) = line(ax, x_years, loc_phox / 1e7, 'color', 'k', 'linewidth', 2, 'marker', 'o', 'markersize', 8);
                l(2) = line(ax, x_years, loc_alpha * 100, 'color', 'k', 'linewidth', 2, 'marker', '^', 'markersize', 8);
                l(3) = line(ax, x_years, loc_vocr, 'color', 'k', 'linewidth', 2, 'marker', 'p', 'markersize', 8);
                
                set(ax, 'xtick', x_years, 'xticklabels', x_labels, 'xlim', [min(x_years)-1, max(x_years)+1], 'xticklabelrotation', 45);
                legend(ax, l, {'P(HO_x) (ppt/s)', '\alpha (%)', 'VOC_R (s^{-1})'});
            end
        end
        
        function plot_epa_nox(cities)
            
            [annual_files, years] = misc_oh_analysis.list_epa_files();
            nox = cell(numel(years), numel(cities));
            for i_yr = 1:numel(annual_files)
                tab = readtable(annual_files(i_yr).name);
                for i_city = 1:numel(cities)
                    xx = misc_oh_analysis.find_epa_sites_for_city(cities{i_city}, tab);
                    pp = strcmpi(tab.ParameterName, 'Oxides of nitrogen (NOx)');
                    
                    nox{i_yr, i_city} = misc_oh_analysis.to_double(tab{xx&pp, {'ArithmeticMean'}});
                end
            end
            
            
            for i_city = 1:numel(cities)
                proto_mat1 = nan(numel(years)-2, 1);
                proto_mat2 = nan(numel(years)-2, 2);
                
                nox_means = proto_mat1;
                nox_medians = proto_mat1;
                nox_quants = proto_mat2;
                for i_yr = 1:size(nox_means,1)
                    yy = i_yr:(i_yr+2);
                    this_nox = veccat(nox{yy, i_city});
                    nox_means(i_yr) = nanmean(this_nox);
                    nox_medians(i_yr) = nanmedian(this_nox);
                    nox_quants(i_yr, :) = quantile(this_nox, [0.05, 0.95]);
                end
                
                l = gobjects(3,1);
                figure;
                [~,l(3)] = plot_error_envelope_y(years(2:end-1), nox_quants(:,1), nox_quants(:,2), [1 0.5 0]);
                l(2) = line(years(2:end-1), nox_medians, 'color', 'r', 'linestyle', '--', 'linewidth', 2);
                l(1) = line(years(2:end-1), nox_means, 'color', 'k', 'linestyle', 'none', 'marker', 'p');
                legend(l, {'Mean', 'Median', '5th/95th percentile'});
                ylabel('[NO_x] (ppb)');
                title(cities{i_city});
            end
            
        end
        
        function plot_epa_voc(cities, varargin)
            p = advInputParser;
            p.addParameter('plot_speciated', true);
            p.parse(varargin{:});
            pout = p.Results;
            plot_speciated = pout.plot_speciated;
            
            [annual_files, years] = misc_oh_analysis.list_epa_files();
            
            proto_mat = nan(numel(years), numel(cities));
            proto_cell = cell(numel(years), numel(cities));
            nmocs = proto_cell;
            vocrs = proto_mat;
            vocr_seoms = proto_mat;
            spec_vocrs = proto_cell;
            
            for i_yr = 1:numel(annual_files)
                tab = readtable(annual_files(i_yr).name);
                for i_city = 1:numel(cities)
                    city_name = cities{i_city};
                    xx = misc_oh_analysis.find_epa_sites_for_city(city_name, tab);
                    pp = strcmp(tab.ParameterName, 'Total NMOC (non-methane organic compound)');
                    this_nmoc = tab{xx&pp, {'ArithmeticMean'}};
                    if iscell(this_nmoc)
                        this_nmoc = cellfun(@str2double, this_nmoc);
                    end
                    nmocs{i_yr, i_city} = this_nmoc;
                    [vocrs(i_yr, i_city), vocr_seoms(i_yr, i_city), spec_vocrs{i_yr, i_city}] = misc_oh_analysis.compute_vocr(tab, city_name);
                end
            end
            
            proto_mat1 = nan(numel(years)-2, numel(cities));
            proto_mat2 = nan(numel(years)-2, numel(cities), 2);
            proto_cell = cell(numel(years)-2, numel(cities));
            nmoc_means = proto_mat1;
            nmoc_errs = proto_mat1;
            nmoc_medians = proto_mat1;
            nmoc_quants = proto_mat2;
            
            vocr_means = proto_mat1;
            vocr_errs = proto_mat1;
            vocr_spec = proto_cell;
            vocr_spec_errs = proto_mat1;
            
            avail_spc = cell(1, numel(cities));
            for i_city = 1:numel(cities)
                avail_spc{i_city} = get_avail_species(spec_vocrs(:, i_city));
            end
            
            for i_yr = 1:size(nmoc_means, 1)
                for i_city = 1:numel(cities)
                    yy = i_yr:(i_yr+2);
                    this_nmoc = veccat(nmocs{yy, i_city});
                    nmoc_means(i_yr, i_city) = nanmean(this_nmoc);
                    nmoc_errs(i_yr, i_city) = nanstd(this_nmoc);
                    nmoc_medians(i_yr, i_city) = nanmedian(this_nmoc);
                    nmoc_quants(i_yr, i_city, :) = quantile(this_nmoc, [0.05, 0.95]);
                    
                    vocr_means(i_yr, i_city) = nanmean(vocrs(yy, i_city));
                    vocr_errs(i_yr, i_city) = sqrt(nanmean(vocr_seoms(yy, i_city).^2));
                    [vocr_spec{i_yr, i_city}, vocr_spec_errs(i_yr, i_city)] = merge_speciated_vocrs(spec_vocrs(yy, i_city), avail_spc{i_city});
                end
            end
            
            for i_city = 1:numel(cities)
                figure;
                title(cities{i_city});
                ax = gca;
                %                 plot_error_envelope_y(years(2:end-1), nmoc_quants(:,1), nmoc_quants(:,2));
                %                 line(years(2:end-1), nmoc_medians, 'color', 'k', 'linewidth', 2);
                %                 line(years(2:end-1), nmoc_means, 'color', 'r', 'marker', 'p');
                %                 ylabel('Total NMOC (ppb C)');
                if plot_speciated
                    spc_mat = cat(1, vocr_spec{:,i_city});
                    total_spc_vocr = sum(spc_mat, 2)';
                    yyaxis(ax, 'left');
                    area(years(2:end-1), spc_mat, 'facecolor', 'flat');
                    colormap(jet);
                    common_opts = {'color', 'k', 'linewidth', 2, 'linestyle', '--'};
                    line(years(2:end-1), total_spc_vocr + vocr_spec_errs(:,i_city)', common_opts{:});
                    line(years(2:end-1), total_spc_vocr - vocr_spec_errs(:,i_city)', common_opts{:});
                    ylabel('VOC_R by species (s^{-1})');
                    
                    yyaxis(ax, 'right');
                    common_opts = {'color', 'r', 'linewidth', 2};
                    line(years(2:end-1), nmoc_means(:,i_city), common_opts{:});
                    line(years(2:end-1), nmoc_means(:,i_city) + nmoc_errs(:,i_city), common_opts{:}, 'linestyle', '--');
                    line(years(2:end-1), nmoc_means(:,i_city) - nmoc_errs(:,i_city), common_opts{:}, 'linestyle', '--');
                    ylabel('Total NMOC (ppbC)');
                else
                    yyaxis(ax, 'left')
                    plot_error_envelope_y(years(2:end-1), vocr_means(:,i_city) - vocr_errs(:,i_city),...
                        vocr_means(:,i_city) + vocr_errs(:,i_city), 'r');
                    
                    yyaxis(ax, 'right')
                    plot_error_envelope_y(years(2:end-1), nmoc_quants(:,1), nmoc_quants(:,2), 'r');
                    
                    yyaxis(ax, 'left')
                    line(years(2:end-1), vocr_means, 'color', 'k', 'linewidth', 2);
                    yylabel('VOC_R (s^{-1})');
                    
                    yyaxis(ax, 'right')
                    
                    line(years(2:end-1), nmoc_medians, 'color', 'r', 'linewidth', 2);
                    line(years(2:end-1), nmoc_means, 'color', 'r', 'marker', 'p', 'linestyle', 'none');
                    ylabel('Total NMOC (ppb C)');
                end
            end
            
            function fns = get_avail_species(spec_vocrs)
                % Since its possible that different years may have
                % different measurements, this will create a struct that
                % has only measurements that are in all time periods
                first = 1;
                while isempty(spec_vocrs{first}.Species)
                    first = first + 1;
                    if first > length(spec_vocrs)
                        fns = {};
                        return
                    end
                end
                fns = spec_vocrs{first}.Species;
                for i=(first+1):numel(spec_vocrs)
                    fns = intersect(fns, spec_vocrs{i}.Species);
                end
            end
            
            function [spc, spc_total_err] = merge_speciated_vocrs(spec_vocrs, fns)
                % Since its possible that different years may have
                % different measurements, this will create a struct that
                % has only measurements that are in all time periods
                
                spc = zeros(1, numel(fns));
                spc_total_err = 0;
                for j=1:numel(fns)
                    spc_vocr = 0;
                    for i=1:numel(spec_vocrs)
                        xx_spc = strcmp(spec_vocrs{i}.Species, fns{j});
                        this_vocr = spec_vocrs{i}.VOCR(xx_spc);
                        this_vocr_seom = spec_vocrs{i}.VOCR_SEOM(xx_spc);
                        if ~isnan(this_vocr)
                            spc_vocr = spc_vocr + this_vocr;
                            spc_total_err = spc_total_err + this_vocr_seom.^2;
                        end
                    end
                    spc(j) = spc_vocr;
                end
                spc_total_err = sqrt(spc_total_err);
            end
            
            
        end
        
        function plot_hcho_refsec_diff(year, dow)
            old_file = misc_emissions_analysis.avg_file_name(year, dow, 'hcho');
            [dirname, basename, ext] = fileparts(old_file);
            new_file = fullfile(dirname, 'NewHCHO', [basename, ext]);
            old_vcds = load(old_file);
            old_vcds = old_vcds.daily;
            new_vcds = load(new_file);
            new_vcds = new_vcds.daily;
            
            tmp = veccat(old_vcds.hcho(:), new_vcds.hcho(:));
            cbrange = calc_plot_limits(tmp(~isoutlier(tmp)));
            tmp = veccat(old_vcds.stddev(:), new_vcds.stddev(:));
            cbrange_sd = calc_plot_limits(tmp(~isoutlier(tmp)));
            
            figure;
            subplot_stretch(2,3);
            
            subplot(2,3,1);
            plot_helper(old_vcds.hcho, cbrange, sprintf('Old HCHO Columns (%d)', year));
            
            subplot(2,3,4);
            plot_helper(old_vcds.stddev, cbrange_sd, 'Std. Dev.');
            
            subplot(2,3,2);
            plot_helper(new_vcds.hcho, cbrange, sprintf('New HCHO Columns (%d)', year));

            subplot(2,3,5);
            plot_helper(new_vcds.hcho, cbrange_sd, 'Std. Dev.');
            
            subplot(2,3,3)
            coldel = new_vcds.hcho - old_vcds.hcho;
            plot_helper(coldel, calc_plot_limits(coldel(:), 'diff'), 'Absolute difference', blue_red_cmap);
            
            subplot(2,3,6);
            colrdel = reldiff(new_vcds.hcho, old_vcds.hcho)*100;
            plot_helper(colrdel, calc_plot_limits(colrdel(~isoutlier(colrdel)), 'diff'), 'Percent difference', blue_red_cmap);
            
            function plot_helper(vcds, cblimits, titlestr, cmap)
                if nargin < 4
                    cmap = parula;
                end
                pcolor(old_vcds.lon, old_vcds.lat, vcds); 
                shading flat
                colorbar;
                colormap(gca, cmap);
                caxis(cblimits);
                title(titlestr);
                state_outlines('k');
            end
        end
        
        function fig = plot_covariance_with_oh(varargin)
            p = advInputParser;
            p.addParameter('loc_inds', 1:71);
            p.addParameter('var1', 'oh');
            p.addParameter('var2', 'vocr'); % any field in the OH type structures
            p.addParameter('corr', true); % true = correlation; false = covariance
            
            p.parse(varargin{:});
            pout = p.Results;
            
            loc_inds = pout.loc_inds;
            main_var = pout.var1;
            co_var = pout.var2;
            do_correlation = pout.corr;
            if do_correlation
                cov_fxn = @corrcoef;
                y_var = 'correlation';
            else
                cov_fxn = @cov;
                y_var = 'covariance';
            end
            
            years = 2006:2013;
            
            oh_types = {'invert', 'invert_hcho', 'invert_hcho_wkday_wkend'};
            type_styles = {'ko', 'b^', 'r*'};
            type_legends = {'Lifetime+SS', 'Lifetime+HCHO SS', 'HCHO SS (VOCR same wkday/wkend)'};
            data = misc_oh_analysis.load_oh_data_for_years(years, 'loc_ids', loc_inds, 'oh_types', oh_types, 'variables', {main_var, co_var});
            var1 = data.(main_var);
            var2 = data.(co_var);
            
            cov_values = nan(numel(loc_inds), numel(oh_types), 2);
            
            for i_loc = 1:numel(loc_inds)
                for i_type = 1:numel(oh_types)
                    for i_dow = 1:2
                        this_oh = var1(:,i_loc,i_type,i_dow);
                        this_var2 = var2(:,i_loc,i_type,i_dow);
                        xx = ~isnan(this_oh) & ~isnan(this_var2);
                        if sum(xx)/numel(this_oh) > 0.5
                            % Only calculate a correlation if we have at
                            % least half of the years
                            tmp_cov = cov_fxn(this_oh(xx), this_var2(xx));
                            cov_values(i_loc, i_type, i_dow) = tmp_cov(1,2); % want the off diagonal term
                        end
                    end 
                end
            end
            
            fig = figure;
            ax1 = subplot(2,1,1);
            hold on
            ax2 = subplot(2,1,2);
            hold on
            
            x = 1:numel(loc_inds);
            for i_type = 1:numel(oh_types)
                plot(ax1, x, cov_values(:, i_type, 1), type_styles{i_type});
                plot(ax2, x, cov_values(:, i_type, 2), type_styles{i_type});
            end
            
            locs = misc_emissions_analysis.read_locs_file();
            loc_inds = misc_emissions_analysis.convert_input_loc_inds(loc_inds);
            locs = misc_emissions_analysis.cutdown_locs_by_index(locs, loc_inds);
            
            legend(ax1, type_legends);
            title(ax1, 'Weekdays');
            legend(ax2, type_legends);
            title(ax2, 'Weekends');
            ax_opts = {'XTickLabel', {locs.ShortName}, 'XTickLabelRotation', 90, 'XLim', [0 x(end)+1], 'YLim', [-1.2, 1.2], 'XTick', x, 'YTick', -1:0.5:1, 'YGrid', 'on'};
            set(ax1, ax_opts{:});
            set(ax2, ax_opts{:});
            
            ylabel_str = sprintf('%s-%s %s', upper(main_var), upper(co_var), y_var);
            ylabel(ax1, ylabel_str);
            ylabel(ax2, ylabel_str);
            subplot_stretch(2,1);
        end
        
        function fig = plot_wrf_hcho_vs_ss_hcho(year, varargin)
            p = advInputParser;
            p.addParameter('phox_mult', 1);
            p.addParameter('alpha_eff', 0.3);
            p.parse(varargin{:});
            pout = p.Results;
            
            alpha_eff = pout.alpha_eff;
            phox_mult = pout.phox_mult;
            
            levels = 1:5;
            ndens = misc_wrf_lifetime_analysis.load_wrf_profiles_for_years(year, 'ndens', 'avg_levels', levels);
            hcho = misc_wrf_lifetime_analysis.load_wrf_profiles_for_years(year, 'hcho', 'avg_levels', levels);
            phox = misc_wrf_lifetime_analysis.load_wrf_profiles_for_years(year, 'phox', 'avg_levels', levels);
            vocr = misc_wrf_lifetime_analysis.load_wrf_profiles_for_years(year, 'vocr', 'avg_levels', levels);
            
            hcho = 1e-6 .* hcho .* ndens;
            phox = 1e-12 .* phox .* ndens;
            
            hcho_ss = hcho_steady_state(vocr, phox * phox_mult, alpha_eff);
            fig = figure;
            scatter(hcho_ss(:), hcho(:), '.');
            plot_fit_line(hcho_ss(:), hcho(:), 'regression', 'rma');
            xlabel('Steady-state HCHO (molec. cm^{-3})')
            ylabel('WRF HCHO (molec. cm^{-3})');
            title(sprintf('\\alpha_{eff} = %.2f, P(HO_x) * %.2f', alpha_eff, phox_mult));
        end
        
        function figs = test_hox_solver_with_wrf(varargin)
            p = advInputParser;
            p.addParameter('loc_inds', 1:71);
            p.addParameter('lifetime_mode', 'no2+oh+ans'); % 'simple', 'complex', 'simple-pre', or 'complex-pre'
            p.addParameter('plot_mode', '1:1'); % 1:1 or time_series
            p.addFlag('three_year_win');
            p.parse(varargin{:});
            pout = p.Results;
            
            % We can test the NO2 + OH -> HNO3 OH calculation and the
            % normal HCHO calculation, but not the weekend/weekday solver
            % because WRF doesn't have a weekend/weekday NOx cycle. 
            %
            % To do so, we need to average OH, NOx, and HCHO concentrations
            % and NOx lifetime from WRF around each location and use the
            % latter three to compute what the steady state model thinks
            % the OH should be. We then compare this against the OH from
            % WRF. 
            %
            % NOx lifetime may be loaded directly or calculated from the
            % loss diagnostics.
            
            three_yr_windows = pout.three_year_win;
            lifetime_mode = pout.lifetime_mode;
            plot_mode = pout.plot_mode;
            
            if three_yr_windows
                years = 2006:2013;
            else
                years = [2005, 2007:2009, 2011:2014];
            end
            n_yrs = numel(years);
            
            loc_inds = misc_emissions_analysis.convert_input_loc_inds(pout.loc_inds);
            locs = misc_emissions_analysis.cutdown_locs_by_index(misc_emissions_analysis.read_locs_file(), loc_inds);
            n_locs = numel(locs);
            
            % Load data for each location we're testing
            hno3_tau_mode = 'hno3';
            switch lower(lifetime_mode)
                case 'simple'
                    lifetime_vars = {'LNOXA','LNOXHNO3'};
                case 'complex'
                    lifetime_vars = {'total_loss','total_prod'};
                case 'simple-pre'
                    lifetime_vars = {'nox_tau_simple'};
                case 'complex-pre'
                    lifetime_vars = {'nox_tau_complex'};
                case 'no2+oh'
                    lifetime_vars = {'no2', 'ho', 'ndens', 'temperature'};
                    hno3_tau_mode = 'no2+oh';
                case 'no2+oh+ans'
                    lifetime_vars = {'LNOXA', 'no2', 'ho', 'ndens', 'temperature'};
                    hno3_tau_mode = 'no2+oh';
                otherwise
                    error('No lifetime variables defined for lifetime_mode = "%s"', lifetime_mode);
            end
            common_vars = {'nox', 'ho', 'LNOXA', 'vocr', 'hcho', 'temperature', 'ndens'};
            calc_vars = {'tau'};
            wrf_vars = unique(veccat(common_vars, lifetime_vars));
            all_vars = veccat(wrf_vars, calc_vars);
            proto_array = repmat({nan(n_yrs,1)}, n_locs, 1);
            % the structure index will be the location, the arrays in each
            % location field will be by year
            wrf_data = make_empty_struct_from_cell(all_vars, proto_array);
            for i_yr=1:n_yrs
                fprintf('Loading data for %d\n', years(i_yr));
                if three_yr_windows
                    year_window = misc_oh_analysis.year_window(years(i_yr));
                    wrf_locs = misc_emissions_analysis.average_profiles_for_locations(year_window, 'TWRF', locs, 'species', wrf_vars);
                    wrf_locs = [wrf_locs.WRFData];
                else
                    year_window = years(i_yr);
                    wrf_locs = misc_oh_analysis.load_wrf_avg(locs, wrf_vars, year_window, []);
                end
                
                for i_loc = 1:n_locs
                    for i_var = 1:numel(wrf_vars)
                        this_var = wrf_vars{i_var};
                        wrf_data(i_loc).(this_var)(i_yr) = wrf_locs(i_loc).(this_var);
                    end
                end
            end
            
            for i_loc = 1:n_locs
                wrf_data(i_loc).tau = misc_oh_analysis.compute_wrf_tau(wrf_data(i_loc), lifetime_mode);
                wrf_data(i_loc).tau_hno3 = misc_oh_analysis.compute_wrf_tau(wrf_data(i_loc), hno3_tau_mode);
                wrf_data(i_loc).tau_ans = misc_oh_analysis.compute_wrf_tau(wrf_data(i_loc), 'ans');
                %wrf_data(i_loc).tau_oh_no2 = misc_oh_analysis.compute_wrf_tau(wrf_data(i_loc), 'no2+oh');
                %(i_loc).tau_oh_no2_ans = misc_oh_analysis.compute_wrf_tau(wrf_data(i_loc), 'no2+oh+ans');
            end
            
            % Calculate the OH solution
            figs = gobjects(n_locs,1);
            for i_loc = 1:n_locs
                figs(i_loc) = figure;
                oh_solutions = make_empty_struct_from_cell({'hno3','invert_hcho'}, proto_array{1});
                vocr_solutions = make_empty_struct_from_cell({'invert_hcho'}, proto_array{1});
                
                
                % Can do the HNO3 solution easily
                k_hno3 = KOHNO2a(wrf_data(i_loc).temperature, wrf_data(i_loc).ndens);
                oh_solutions.hno3 = 1 ./ (k_hno3 .* wrf_data(i_loc).tau .* 3600);
                this_nox = wrf_data(i_loc).nox .* wrf_data(i_loc).ndens .* 1e-6;
                this_hcho = wrf_data(i_loc).hcho .* wrf_data(i_loc).ndens .* 1e-6;
                
                oh_hcho = oh_solutions.invert_hcho;
                vocr_hcho = vocr_solutions.invert_hcho;
                phox_ss = proto_array{1};
                tau_ss = proto_array{1};
                tau_hno3_ss = proto_array{1};
                tau_ans_ss = proto_array{1};
                parfor i_yr = 1:n_yrs
                    fprintf('Solving OH and VOCR for %d\n', years(i_yr));
                    [oh, ~, ~, soln] = hox_solve_tau_hcho_constraint(this_nox(i_yr), wrf_data(i_loc).tau(i_yr), this_hcho(i_yr), 0.04);
                    oh_hcho(i_yr) = oh;
                    vocr_hcho(i_yr) = soln.vocr;
                    phox_ss(i_yr) = soln.phox ./ (2e19*1e-12); %convert ~ to ppt/s to put on ~ same scale as vocr
                    tau_ss(i_yr) = soln.tau;
                    tau_hno3_ss(i_yr) = soln.tau_hno3;
                    tau_ans_ss(i_yr) = soln.tau_ans;
                end
                oh_solutions.invert_hcho = oh_hcho;
                vocr_solutions.invert_hcho = vocr_hcho;
                
                % there's a bit of extra back-and-forth into the
                % structure/out of the structure because of old refactoring
                oh_ss = oh_solutions.invert_hcho;
                oh_hno3 = oh_solutions.hno3;
                oh_wrf = wrf_data(i_loc).ho .* wrf_data(i_loc).ndens .* 1e-6;
                vocr_ss = vocr_solutions.invert_hcho;
                vocr_wrf = wrf_data(i_loc).vocr;
                tau_wrf = wrf_data(i_loc).tau;
                tau_hno3_wrf = wrf_data(i_loc).tau_hno3;
                tau_ans_wrf = wrf_data(i_loc).tau_ans;
                
                if strcmpi(plot_mode, 'time_series')
                    common_opts = {'markersize', 12, 'linewidth', 2};
                    subplot(2,1,1);
                    yyaxis left;
                    l = gobjects(6,1);
                    l(1) = line(years, oh_wrf, 'color', 'k', 'marker', '^', common_opts{:});
                    l(2) = line(years, oh_hno3, 'color', 'k', 'marker', 'o', common_opts{:});
                    l(3) = line(years, oh_hcho, 'color', 'k', 'marker', '*', common_opts{:});
                    ylabel('[OH] (molec. cm^{-3}');
                    
                    
                    yyaxis right;
                    l(4)=line(years, vocr_wrf, 'color', 'r', 'marker', '^', 'linestyle', 'none', common_opts{:});
                    l(5)=line(years, vocr_ss, 'color', 'r', 'marker', '*', 'linestyle', 'none', common_opts{:});
                    l(6)=line(years, phox_ss, 'color', 'b', 'marker', '*', 'linestyle', 'none', common_opts{:});
                    ylabel('VOC_R (s^{-1}) and P(HO_x) (ppt s^{-1})')
                    
                    legend(l, {'WRF OH', 'HNO_3 OH', 'HCHO SS OH', 'WRF VOC_R', 'HCHO SS VOC_R', 'HCHO SS P(HO_x)'}, 'location', 'eastoutside');
                    title(locs(i_loc).ShortName);
                    
                    subplot(2,1,2);
                    l2 = gobjects(6,1);
                    l2(1) = line(years, tau_wrf, 'color', 'k', 'marker', 'o', common_opts{:});
                    l2(2) = line(years, tau_hno3_wrf, 'color', 'k', 'marker', '^', common_opts{:});
                    l2(3) = line(years, tau_ans_wrf, 'color', 'k', 'marker', '*', common_opts{:});
                    l2(4) = line(years, tau_ss, 'color', 'r', 'marker', 'o', common_opts{:});
                    l2(5) = line(years, tau_hno3_ss, 'color', 'r', 'marker', '^', common_opts{:});
                    l2(6) = line(years, tau_ans_ss, 'color', 'r', 'marker', '*', common_opts{:});
                    legend(l2, {'WRF total', 'WRF HNO_3', 'WRF ANs', 'SS total', 'SS HNO_3', 'SS ANs'}, 'location', 'eastoutside');
                    ylabel('\tau (h)');
                    
                    subplot_stretch(2,1.5)
                elseif strcmpi(plot_mode, '1:1')
                    subplot(2,1,1);
                    hold on
                    l1 = gobjects(4,1);
                    l1(1) = scatter(oh_wrf, oh_hno3, [], years, '^');
                    l1(3) = scatter(oh_wrf, oh_ss, [], years, 'o', 'filled');
                    [x,y,lstr_hno3] = calc_fit_line(oh_wrf, oh_hno3, 'regression', 'rma');
                    l1(2) = line(x,y,'color', 'k', 'linestyle', '--', 'linewidth', 2);
                    [x,y,lstr_hcho] = calc_fit_line(oh_wrf, oh_ss, 'regression', 'rma');
                    l1(4) = line(x,y,'color', 'b', 'linestyle', '-.', 'linewidth', 2);
                    legend(l1, {'HNO3', lstr_hno3, 'HCHO SS', lstr_hcho});
                    cb = colorbar;
                    cb.Label.String = 'Year';
                    colormap copper
                    xlabel('WRF OH (molec. cm^{-3})');
                    ylabel('Solved OH (molec. cm^{-3})');
                    title(locs(i_loc).ShortName);
                    
                    subplot(2,1,2);
                    hold on
                    l2 = gobjects(2,1);
                    l2(1) = scatter(vocr_wrf, vocr_ss, [], years, 'filled');
                    [x,y,lstr_vocr] = calc_fit_line(vocr_wrf, vocr_ss, 'regression', 'rma');
                    l2(2) = line(x, y, 'color', 'r', 'linestyle', '--', 'linewidth', 2);
                    legend(l2(2), {lstr_vocr});
                    cb = colorbar;
                    cb.Label.String = 'Year';
                    colormap copper
                    xlabel('WRF VOC_R (s^{-1})');
                    ylabel('Solved VOC_R (s^{-1})');
                    subplot_stretch(2,1);
                end
            end
            
            
        end
        
        function tau = compute_wrf_tau(Wrf, lifetime_mode)
            switch lower(lifetime_mode)
                case 'simple'
                    loss = Wrf.LNOXHNO3 + Wrf.LNOXA;
                    nox = get_nox(Wrf);
                    tau = nox ./ loss ./ 3600;
                case 'ans'
                    loss = Wrf.LNOXA;
                    nox = get_nox(Wrf);
                    tau = nox ./ loss ./ 3600;
                case 'hno3'
                    loss = Wrf.LNOXHNO3;
                    nox = get_nox(Wrf);
                    tau = nox ./ loss ./ 3600;
                case 'no2+oh'
                    loss = misc_oh_analysis.get_loss_via_oh_no2(Wrf);
                    nox = get_nox(Wrf);
                    tau = nox ./ loss ./ 3600;
                case 'no2+oh+ans'
                    loss = misc_oh_analysis.get_loss_via_oh_no2(Wrf) + Wrf.LNOXA;
                    nox = get_nox(Wrf);
                    tau = nox ./ loss ./ 3600;
                case 'complex'
                    loss = Wrf.total_loss - Wrf.total_prod;
                    nox = get_nox();
                    tau = nox ./ loss ./ 3600;
                case 'simple-pre'
                    tau = Wrf.nox_tau_simple;
                case 'complex-pre'
                    tau = Wrf.nox_tau_complex;
                otherwise
                    error('No tau computation for lifetime_mode = "%s"', lifetime_mode);
            end
            
            function nox = get_nox(Wrf)
                if isfield(Wrf, 'nox')
                    nox = Wrf.nox;
                else
                    nox = Wrf.no + Wrf.no2;
                end
            end
            
            
        end
        
        function value = get_loss_via_oh_no2(Wrf)
            T = Wrf.temperature;
            ndens = Wrf.ndens;
            k = wrf_rate_expr('TROE', 1.49e-30 , 1.8 , 2.58e-11 , 0.0, T, ndens);
            oh_nd = Wrf.ho .* 1e-6 .* ndens;
            no2_nd = Wrf.no2 .* 1e-6 .* ndens;
            loss_nd = k .* oh_nd .* no2_nd;
            value = loss_nd ./ ndens .* 1e6; % return in ppm/s
        end
        
        function compare_hox_ss_models(varargin)
            p = advInputParser;
            p.addParameter('var1','nox'); % this and var2 can be nox, phox, vocr, or alpha
            p.addParameter('var2','phox');
            p.parse(varargin{:});
            pout = p.Results;
            
            var1 = pout.var1;
            var2 = pout.var2;
            
            ppbMolec = 2.45e10;
            defaults.nox = logspace(-1,1,20)*2; % in ppb
            defaults.vocr = 1:20;
            defaults.phox = 5e6:1e6:5e7;
            defaults.alpha = 0.01:0.01:0.1;
            
            fns = fieldnames(defaults);
            xx = ~ismember(fns, {var1, var2});
            var34 = fns(xx);
            
            [x,y] = meshgrid(defaults.(var1), defaults.(var2));
            z1 = repmat(mean(defaults.(var34{1})), size(x));
            z2 = repmat(mean(defaults.(var34{2})), size(x));
            vars = veccat({var1, var2}, var34, 'column');
            values = {x, y, z1, z2};
            
            nox = values{strcmp(vars,'nox')};
            vocr = values{strcmp(vars,'vocr')};
            phox = values{strcmp(vars, 'phox')};
            alpha = values{strcmp(vars, 'alpha')};
            
            oh_a = nan(size(x));
            oh_b = nan(size(x));
            parfor i=1:numel(x)
                fprintf('OH %d of %d\n', i, numel(x));
                oh_a(i) = hox_ss_solver(nox(i) * ppbMolec, phox(i), vocr(i), alpha(i));
                [~, ~, oh_b(i)] = Analytic_P_O3_wNOx(nox(i), phox(i), vocr(i), alpha(i), 4);
            end
            
            figure;
            oh_diff = oh_a - oh_b;
            contourf(x,y,oh_diff);
            cb = colorbar;
            caxis(calc_plot_limits(oh_diff(:),'diff'));
            colormap(blue_red_cmap);
            cb.Label.String = '\Delta [OH] (molec. cm^{-3}, solver - Murphy)';
            xlabel(upper(var1));
            ylabel(upper(var2));
        end
        
        function [wrf_locs, figs] = test_hox_ss_with_wrf(varargin)
            % What we want to do is get the NOx, VOCR, P(HOx) and alpha
            % from WRF, plug these into HOX_SS_SOLVER and see how the OH it
            % produces compares with the OH from WRF itself.
            %
            % We want the option to use either time-average or
            % instantaneous WRF files, and to do the averaging within
            % different radii around the locations specified.
            
            p = advInputParser;
            p.addParameter('loc_inds', 1:71);
            p.addParameter('avg_radius', []); % radius around cities to include. An empty array uses box width, a number is in degrees
            p.addParameter('time_mode', 'avg'); % 'avg' or 'inst'
            p.addParameter('local_hour', 13); % local hour to use for 'inst' mode
            p.addParameter('years', [2005, 2007:2009, 2011:2014]); % what years to use. if 'inst', see next option. Skip 2006 and 2010 by default b/c they lack full WRF output
            p.addParameter('inst_dates', {'04-01','09-30'}); % what days of year (in mm-dd format) to use for 'inst'
            p.addParameter('oh_model','solver'); % 'solver' or 'murphy'
            p.addParameter('vocr_mult', 1);
            
            p.parse(varargin{:});
            pout = p.Results;
            
            loc_inds = misc_emissions_analysis.convert_input_loc_inds(pout.loc_inds);
            avg_radius = pout.avg_radius;
            time_mode = pout.time_mode;
            local_hour = pout.local_hour;
            years = pout.years;
            inst_dates = pout.inst_dates;
            oh_model = pout.oh_model;
            vocr_mult = pout.vocr_mult;
            
            % First need to get the locations
            locs = misc_emissions_analysis.cutdown_locs_by_index(misc_emissions_analysis.read_locs_file(), loc_inds);
            
            % Next load the data, averaging around each location
            wrf_vars = {'no', 'no2', 'ho', 'vocr', 'alpha', 'phox', 'ndens', 'LNOXA', 'LNOXHNO3'};
            switch lower(time_mode)
                case 'avg'
                    wrf_locs = misc_oh_analysis.load_wrf_avg(locs, wrf_vars, years, avg_radius);
                case 'inst'
                    wrf_locs = load_wrf_inst(locs);
                otherwise
                    error('time_mode "%s" not recognized', time_mode);
            end
            
            % Compute what the steady state solver thinks the OH should be.
            % Plot how this compares to the WRF OH
            figs = gobjects(numel(locs),1);
            if isempty(avg_radius)
                radius_str = 'by location';
            else
                radius_str = sprintf('%.2f deg', avg_radius);
            end
            for j_loc = 1:numel(locs)
                switch lower(time_mode)
                    case 'avg'
                        col_var = years;
                        col_label = 'Year';
                    case 'inst'
                        col_var = wrf_locs(j_loc).PHOTR_NO2 ./ 60;
                        col_label = 'j(NO_2) (s^{-1})';
                    otherwise
                        error('No color variable defined for time mode = "%s"', time_mode)
                end
                wrf_locs(j_loc).ss_oh = solve_oh(wrf_locs(j_loc));
                wrf_locs(j_loc).wrf_tau = misc_oh_analysis.compute_wrf_tau(wrf_locs(j_loc), 'simple');
                wrf_locs(j_loc).ss_tau = solve_tau(wrf_locs(j_loc));
                
                figs(j_loc) = figure;
                subplot(2,1,1);
                oh_wrf = wrf_locs(j_loc).ho .* wrf_locs(j_loc).ndens * 1e-6;
                oh_ss = wrf_locs(j_loc).ss_oh .* wrf_locs(j_loc).ndens * 1e-6;
                scatter(oh_wrf, oh_ss, [], col_var);
                plot_fit_line(oh_wrf, oh_ss, 'one2one', false, 'regression', 'rma', 'legend_loc', 'best');
                cb = colorbar;
                colormap copper
                cb.Label.String = col_label;
                xlabel('WRF [OH] (molec cm^{-3})');
                ylabel('Steady state [OH] (molec cm^{-3})');
                title(sprintf('%s (radius %s)', locs(j_loc).ShortName, radius_str));
                
                subplot(2,1,2);
                scatter(wrf_locs(j_loc).wrf_tau, wrf_locs(j_loc).ss_tau, [], col_var);
                plot_fit_line(wrf_locs(j_loc).wrf_tau, wrf_locs(j_loc).ss_tau, 'one2one', false, 'regression', 'rma', 'legend_loc', 'best');
                cb = colorbar;
                colormap copper
                cb.Label.String = col_label;
                xlabel('WRF NO_x lifetime (h)');
                ylabel('Steady state NOx lifetime (h)');
            end
            
            wrf_locs = copy_structure_fields(locs, wrf_locs, 'missing');
            
            function wrf_data = load_wrf_inst(locs)
                % 1) Figure out what hours need loaded
                % 2) Create date vector and loop over days
                % 3) Load the variables from the required hours and average
                %    the right hour for the right city. PHOx, VOCR, and
                %    alpha will all need calculated.
                
                % 1) figure out hours
                loc_long = [locs.Longitude];
                local_hrs = local_hour - round(loc_long / 15);
                local_hrs_to_load = unique(local_hrs);
                
                % 2) start looping over days
                sdates = arrayfun(@(y) sprintf('%04d-%s', y, inst_dates{1}), years, 'uniform', false);
                edates = arrayfun(@(y) sprintf('%04d-%s', y, inst_dates{2}), years, 'uniform', false);
                dvec = make_datevec(sdates, edates);
                
                proto_array = repmat({nan(numel(dvec),1)},size(locs));
                wrf_data = make_empty_struct_from_cell(veccat({'PHOTR_NO2'}, wrf_vars), proto_array);
                
                last_year = nan;
                for i_day = 1:numel(dvec)
                    fprintf('Loading %.2f%%\n', (i_day-1)./numel(dvec)*100);
                    this_year = year(dvec(i_day));
                    if last_year ~= this_year
                        % get the right processing function for each year;
                        % since some years can't calculate these quantities
                        last_year = this_year;
                        alpha_proc = misc_wrf_lifetime_analysis.setup_alpha_calc(this_year);
                        phox_proc = misc_wrf_lifetime_analysis.setup_phox_calc(this_year);
                        vocr_proc = misc_wrf_lifetime_analysis.setup_vocr_calc(this_year);
                        % to avoid double-loading certain variables (e.g. ndens)
                        % make a list of all unique variables to load from each
                        % file
                        vars_to_load = unique(veccat({'PHOTR_NO2'}, wrf_vars, alpha_proc.variables, phox_proc.variables, vocr_proc.variables, 'column'));
                        xx = ~ismember(vars_to_load, {'alpha','phox','vocr'});
                        vars_to_load = vars_to_load(xx);
                    end
                    
                    wrf_files = cell(size(local_hrs_to_load));
                    for i_hr = 1:numel(wrf_files)
                        wrf_datetime = dvec(i_day) + local_hrs_to_load(i_hr)/24;
                        wrf_files{i_hr} = find_wrf_path('us','daily',wrf_datetime,'fullpath');
                    end
                    
                    if i_day == 1
                        % the coordinates are the same across all files
                        xlon = ncread(wrf_files{1}, 'XLONG');
                        xlat = ncread(wrf_files{1}, 'XLAT');
                    end
                    
                    wrf_profiles = read_wrf_vars('', wrf_files, vars_to_load, 'squeeze', 1, 'as_struct');
                   
                    wrf_vars = veccat({'PHOTR_NO2'}, wrf_vars, 'column');
                    for i_var = 1:numel(wrf_vars)
                        this_var = wrf_vars{i_var};
                        switch lower(this_var)
                            case 'alpha'
                                data = alpha_proc.proc_fxn(wrf_profiles);
                            case 'phox'
                                data = phox_proc.proc_fxn(wrf_profiles);
                            case 'vocr'
                                data = vocr_proc.proc_fxn(wrf_profiles);
                            otherwise
                                % Most of the variables can be read
                                % by read_wrf_vars
                                data = wrf_profiles.(this_var);
                        end
                        
                        % Now we need to average the data for each location
                        for i_loc = 1:numel(locs)
                            % data will be nlon x nlat x nlevels x nhours
                            % so get which slice in the last dimension we
                            % should use
                            i_hr = local_hrs(i_loc) == local_hrs_to_load;
                            wrf_data(i_loc).(this_var)(i_day) = misc_wrf_lifetime_analysis.average_wrf_data_around_loc(locs(i_loc),...
                                data(:,:,:,i_hr), xlon, xlat, 'avg_levels', 1:5, 'radius', avg_radius);
                        end
                    end
                end
                
                fprintf('Loading 100%%\n');
            end
            
            
            function oh = solve_oh(wrf_data)
                switch lower(oh_model)
                    case 'solver'
                        nox = (wrf_data.no + wrf_data.no2) .* 1e-6 .* wrf_data.ndens;
                        oh_calc_fxn = @hox_ss_solver;
                    case 'murphy'
                        nox = (wrf_data.no + wrf_data.no2) .* 1e3;
                        oh_calc_fxn = @Analytic_P_O3_wNOx;
                    otherwise
                        error('oh_model "%s" not recognized');
                end
                ratio = wrf_data.no2 ./ wrf_data.no;
                vocr = wrf_data.vocr * vocr_mult;
                % phox in ppt/s -> molec/cm3/s
                phox = wrf_data.phox .* wrf_data.ndens .* 1e-12;
                alpha = 2 ./ (1./wrf_data.alpha+2); % calculated alpha as LNOXA./PO3, which is slightly wrong
                oh = nan(size(nox));
                for i = 1:numel(oh)
                    try
                        oh(i) = oh_calc_fxn(nox(i), phox(i), vocr(i), alpha(i), ratio(i));
                    catch err
                        if strcmp(err.identifier, 'hox_ss_solver:invalid_input')
                            continue  % input invalid (probably a NaN). Can't compute OH
                        else
                            rethrow(err)
                        end
                    end
                end
                oh = oh ./ wrf_data.ndens .* 1e6; % return in ppm just like WRF
            end
            
            function tau = solve_tau(wrf_data)
                nox = (wrf_data.no + wrf_data.no2) .* 1e-6 .* wrf_data.ndens;
                no2_no = wrf_data.no2 ./ wrf_data.no;
                vocr = wrf_data.vocr * vocr_mult;
                % phox in ppt/s -> molec/cm3/s
                phox = wrf_data.phox .* wrf_data.ndens .* 1e-12;
                alpha = 2 ./ (1./wrf_data.alpha+2); % calculated alpha as LNOXA./PO3, which is slightly wrong
                tau = nan(size(nox));
                for i = 1:numel(tau)
                    tau(i) = nox_lifetime(nox(i), 'phox', phox(i), 'alpha', alpha(i), 'vocr', vocr(i), 'no2_no', no2_no(i));
                end
            end
        end
        %%%%%%%%%%%%%%%%%%%%%%%%%
        % Plot helper functions %
        %%%%%%%%%%%%%%%%%%%%%%%%%
        
        function col = get_time_period_color(time_period)
            time_period = mean(time_period);
            all_tps = cellfun(@mean, misc_oh_analysis.all_time_periods);
            col = map2colmap(time_period, min(all_tps), max(all_tps), jet);
        end
        
        function [vocr, vocr_sigma, speciated_vocr] = compute_vocr(epa_table, city)
            % COMPUTE_VOCR Calculate the VOCR from EPA insitu measurements
            %
            % [VOCR, VOCR_SIGMA] = COMPUTE_VOCR(EPA_TABLE, CITY) computes
            % the VOCR and its uncertainty from the table EPA_TABLE
            % (created by loading one of the
            % "annual_conc_by_monitor_yyyy.csv" files with READTABLE) for
            % the city CITY.
            %
            % VOCR is the sum of the average VOCRs for each measured
            % species averaged across all sites in CITY. VOCR_SIGMA is the
            % quadrature sum of the standard errors of the mean for each
            % species across all measurement sites.
            
            % For each species listed in the database, find it in the
            % table, then figure out if its in units of ppbC or ng/m^3. The
            % former is defined by
            % https://www3.epa.gov/ttnamti1/files/ambient/npapsop/sop017.pdf
            % as ppbV * number of carbons (pg. 11). So to get ppb on the
            % way to number density, we need to divide by the number of
            % carbons. We'll handle that by counting carbons in the SMILES
            % string. For the others, we just need the molecular weight to
            % convert. But there's no point in doing any of this if we
            % didn't find an MCM rate constant, so check that first.
            
            db = misc_oh_analysis.rate_db;
            tmp_data = db.fetch('select Parameter, SMILES, MCMRate, MolecularWeight from rates');
            species = tmp_data(:,1);
            smiles = tmp_data(:,2);
            rates = tmp_data(:,3);
            mol_wt = tmp_data(:,4);
            
            vocr = 0;
            vocr_sigma = 0;
            speciated_vocr = cell(0,3);
            missing_species = {};
            
            ll = misc_oh_analysis.find_epa_sites_for_city(city, epa_table);
            
            for i_sp = 1:numel(species)
                if length(rates{i_sp}) <= 1
                    % NULL values for mising rates show up as empty
                    % strings; also some of them I had marked with a "!" to
                    % go back and double check, this skips both.
                    continue
                end
                
                oh_rxn_rate = eval_rate(rates{i_sp});
                
                % Now get the concentrations and units for the current 
                % parameter
                xx = strcmp(epa_table.ParameterName, species{i_sp}) & ll;
                if sum(xx) == 0
                    missing_species{end+1} = species{i_sp};
                    continue
                end
                
                % Double check that we didn't accidentally get sites in two
                % very different cities. Use 2 degree as the cutoff b/c
                % that's how wide the EMG boxes are.
                max_separation = calc_max_dist(misc_oh_analysis.to_double(epa_table.Longitude(xx)),...
                    misc_oh_analysis.to_double(epa_table.Latitude(xx)));
                if max_separation > 2
                    error('Sites chosen for "%s" have a maximum separation of %f degrees - exceeds the maximum allowed', city, max_separation)
                end
                
                species_concs = misc_oh_analysis.to_double(epa_table.ArithmeticMean(xx));
                species_units = unique(epa_table.UnitsOfMeasure(xx));
                if numel(species_units) ~= 1
                    error('Multiple units found (%s) for species "%s" in city "%s"', strjoin(species_units, ', '), species{i_sp}, city);
                else
                    species_units = species_units{1};
                end
                
                switch species_units
                    case 'Parts per billion Carbon'
                        n_carbon = count_carbons(smiles{i_sp});
                        % convert ppbC -> ppb -> molec/cm^3
                        species_concs = species_concs ./ n_carbon .* 1e-9 .* 2e19;
                    case 'Nanograms/cubic meter (25 C)'
                        mw = mol_wt{i_sp};
                        % ng / m^3 -> g / m^3 -> mol/m^3 -> molec./m^3 ->
                        % molec/cm^3
                        species_concs = species_concs .* 1e-9 ./ mw .* 6.022e23 ./ (100^3);
                    otherwise
                        error('No conversion defined for units "%s"', species_units);
                end
                
                % Finally VOCR is the concentration of the VOC in number
                % density times its OH rate constant. Since some species
                % may not be measured at all sites, we'll average each
                % species and compute its standard error of the mean. We'll
                % add the errors in quadrature to get the total uncertainty
                % in the VOCR. Using SEOM to accound for the fact that more
                % sites measuring a particular species reduce the
                % uncertainty of that species.
                
                this_vocr = species_concs .* oh_rxn_rate;
                this_vocr_seom = nanstd(this_vocr) ./ sum(~isnan(this_vocr));
                this_vocr = nanmean(this_vocr);
                
                if ~isnan(this_vocr)  % just in case all concentrations were NaN
                    vocr = vocr + this_vocr;
                    vocr_sigma = vocr_sigma + this_vocr_seom.^2;
                    speciated_vocr = cat(1, speciated_vocr, {species{i_sp}, this_vocr, this_vocr_seom});
                end
            end
            
            if isempty(speciated_vocr)
                % If nothing was found, set VOCR to NaN, not 0.
                vocr = nan;
                vocr_sigma = nan;
            end
            
            vocr_sigma = sqrt(vocr_sigma);
            speciated_vocr = cell2table(speciated_vocr, 'VariableNames', {'Species', 'VOCR', 'VOCR_SEOM'});
            if numel(missing_species) > 0
                missing_species_list = strjoin(missing_species, '\n  * ');
                warning('%d of %d species could not be found for %s:\n  * %s', numel(missing_species), numel(species), city, missing_species_list);
            end
            
            function nc = count_carbons(smiles)
                % Look for a C not followed by a lower case letter (e.g. Cl
                % = chlorine, not a carbon). Include the end of the
                % string.
                nc = numel(regexp(smiles, 'C(?=[^a-z])|C$'));
            end
            
            function rate = eval_rate(rate_expr)
                rate_expr = strrep(rate_expr, '@', '^');
                rate_expr = strrep(rate_expr, 'EXP', 'exp');
                TEMP = 298;
                M = 2e19;
                
                switch rate_expr
                    case 'KMT15'
                        rate = KMT15(TEMP, M);
                    case 'KMT16'
                        rate = KMT16(TEMP, M);
                    case 'KMT17'
                        rate = KMT17(TEMP, M);
                    otherwise
                        rate = eval(rate_expr);
                end
            end
            
            function r_longest = calc_max_dist(lon, lat)
                r_mat = zeros(numel(lon));
                for i=1:numel(lon)
                    for j=(i+1):numel(lon)
                        r = (lon(i) - lon(j)).^2 + (lat(i) - lat(j)).^2;
                        r_mat(i,j) = r;
                    end
                end
                r_longest = max(sqrt(r_mat(:)));
            end
        end
        
        function data = load_oh_data_for_years(years, varargin)
            p = advInputParser;
            p.addParameter('loc_ids', 1:71);
            p.addParameter('oh_types', {}); % leave empty for all
            p.addParameter('variables', {'oh'}); % cell array of fieldnames in each OH struct. If missing, then will be skipped
            
            p.parse(varargin{:});
            pout = p.Results;
            
            loc_inds = misc_emissions_analysis.convert_input_loc_inds(pout.loc_ids);
            oh_types = pout.oh_types;
            variables = pout.variables;
            
            for i_yr = 1:numel(years)
                yr = years(i_yr);
                yr_window = (yr-1):(yr+1);
                fprintf('Loading %s\n', sprintf_ranges(yr_window));
                OH = load(misc_emissions_analysis.oh_file_name(yr_window));
                OH.locs_wkday = misc_emissions_analysis.cutdown_locs_by_index(OH.locs_wkday, loc_inds);
                OH.locs_wkend = misc_emissions_analysis.cutdown_locs_by_index(OH.locs_wkend, loc_inds);
                
                if i_yr == 1
                    if isempty(oh_types)
                        oh_types = fieldnames(OH.locs_wkday(1).OH);
                    end
                    % each array will be years x locs x oh type x days of
                    % week
                    data = make_empty_struct_from_cell(variables, nan(numel(years), numel(loc_inds), numel(oh_types), 2));
                end
                
                for i_loc = 1:numel(loc_inds)
                    for i_type = 1:numel(oh_types)
                        this_type = oh_types{i_type};
                        weekday_oh = OH.locs_wkday(i_loc).OH.(this_type);
                        weekend_oh = OH.locs_wkend(i_loc).OH.(this_type);
                        for i_var = 1:numel(variables)
                            this_var = variables{i_var};
                            if isfield(weekday_oh, this_var)
                                data.(this_var)(i_yr, i_loc, i_type, 1) = weekday_oh.(this_var);
                            end
                            if isfield(weekend_oh, this_var)
                                data.(this_var)(i_yr, i_loc, i_type, 2) = weekend_oh.(this_var);
                            end
                        end
                    end
                end
            end
        end
        
        function wrf_data = load_wrf_avg(locs, wrf_vars, years, avg_radius)
            avg_levels = 1:5;
            proto_array = repmat({nan(numel(years),1)},size(locs));
            wrf_data = make_empty_struct_from_cell(wrf_vars, proto_array);
            for i_loc = 1:numel(locs)
                fprintf('Loading: %.1f%%\n', (i_loc - 1)/numel(locs)*100);
                for i_yr = 1:numel(years)
                    for i_var = 1:numel(wrf_vars)
                        this_var = wrf_vars{i_var};
                        if strcmpi(this_var, 'nox')
                            no = misc_wrf_lifetime_analysis.average_profiles_around_loc(locs(i_loc),...
                                years(i_yr), 'no', 'radius', avg_radius, 'avg_levels', avg_levels);
                            no2 = misc_wrf_lifetime_analysis.average_profiles_around_loc(locs(i_loc),...
                                years(i_yr), 'no2', 'radius', avg_radius, 'avg_levels', avg_levels);
                            wrf_data(i_loc).(this_var)(i_yr) = no + no2;
                        else
                            wrf_data(i_loc).(this_var)(i_yr) = misc_wrf_lifetime_analysis.average_profiles_around_loc(locs(i_loc),...
                                years(i_yr), this_var, 'radius', avg_radius, 'avg_levels', avg_levels);
                        end
                    end
                end
            end
        end
        
        %%%%%%%%%%%%%%%%
        % Misc Methods %
        %%%%%%%%%%%%%%%%
        
        function win = year_window(yr)
            win = (yr-1):(yr+1);
        end
        
        function rows = get_rates()
            fid = fopen(fullfile(misc_oh_analysis.epa_data_dir, 'rate4.csv'));
            tline = fgetl(fid);
            fields = strsplit(tline, '\t');
            tline = fgetl(fid);
            rows = struct([]);
            while ischar(tline)
                tline = regexprep(tline, '\%.+$', '');
                if isempty(tline)
                    tline = fgetl(fid);
                    continue
                end
                line_data = strsplit(tline, '\t');
                line_data = cat(2, line_data, repmat({''}, 1, numel(fields)-numel(line_data)));
                this_line = make_struct_from_field_values(fields, line_data);
                if isempty(rows)
                    rows = this_line;
                else
                    rows(end+1) = this_line;
                end
                tline = fgetl(fid);
            end
            fclose(fid);
        end
        
        function result = get_param_info()
            mcm = misc_oh_analysis.get_rates();
            tab = readtable(fullfile(misc_oh_analysis.epa_data_dir, 'annual_conc_by_monitor_2005.csv'));
            
            params = {'MethodName', 'MetricUsed', 'UnitsOfMeasure'};
            values = struct([]);
            
            for i_var = 1:numel(mcm)
                this_param = mcm(i_var).Parameter;
                xx = strcmp(tab.ParameterName, this_param);
                if sum(xx) == 0
                    fprintf('Could not find %s in the original table\n', this_param);
                    continue
                end
                
                values(end+1).Parameter = this_param;
                
                for p = 1:numel(params)
                    rows = unique(tab{xx, params(p)});
                    if numel(rows) > 1
                        strval = sprintf('(%s)', strjoin(rows, ', '));
                    else
                        strval = rows{1};
                    end
                    values(end).(params{p}) = strval;
                end
            end
            result = struct2table(values);
        end
        
        function xx = find_epa_sites_for_city(city, epa_table)
            loc = misc_oh_analysis.match_loc_to_city(city);
            xx = misc_oh_analysis.find_epa_sites_for_loc(loc, epa_table);
        end
        
        function xx = find_epa_sites_for_loc(loc, epa_table)
            lon = loc.Longitude;
            lat = loc.Latitude;
            radius = mean(loc.BoxSize(3:4));
            
            site_lons = misc_oh_analysis.to_double(epa_table.Longitude);
            site_lats = misc_oh_analysis.to_double(epa_table.Latitude);
            
            xx = (site_lons - lon).^2 + (site_lats - lat).^2 < radius .^2;
        end
        
        function val = to_double(val)
            if iscell(val)
                val = cellfun(@str2double, val);
            end
        end
        
        function loc = match_loc_to_city(city)
            all_locs = misc_emissions_analysis.read_locs_file();
            for i = 1:numel(all_locs)
                if strcmpi(all_locs(i).Location, city) || strcmpi(all_locs(i).ShortName, city)
                    loc = all_locs(i);
                    return
                end
            end
            
            error('Could not find a location matching "%s"', city);
        end
    end
end

