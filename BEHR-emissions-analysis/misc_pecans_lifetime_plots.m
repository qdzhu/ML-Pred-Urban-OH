classdef misc_pecans_lifetime_plots
    %UNTITLED Summary of this class goes here
    %   Detailed explanation goes here
    
    properties
        Property1
    end
    
    methods(Static)
        function value = workspace_dir()
            my_dir = fileparts(mfilename('fullpath'));
            value = fullfile(my_dir, 'Workspaces');
        end
        
        function value = pecans_dir()
            value = fullfile(misc_pecans_lifetime_plots.workspace_dir, 'PECANS');
        end
        
        
        
        function plot_normalized_comparison
            % Compared the normalized and centered line densities from a
            % city's BEHR data and the PECANS outputs
            
            pecans_data = misc_pecans_lifetime_plots.read_pecans_for_tau(3);
            behr_data = load(misc_emissions_analysis.fits_file_name('2012-04-01','2014-09-30',false,1:71,'UMTWRFS','lu'));
            behr_loc = behr_data.locs(9); % use Chicago
            
            % Do need to transform PECANS x from meters to kilometers
            [pecans_x, pecans_ld] = misc_pecans_lifetime_plots.center_and_normalize(pecans_data.x, pecans_data.A);
            [behr_x, behr_ld] = misc_pecans_lifetime_plots.center_and_normalize(behr_loc.no2_sectors.x, behr_loc.no2_sectors.linedens);
            
            figure; plot(behr_x, behr_ld, 'b', pecans_x, pecans_ld, 'r');
            legend(sprintf('BEHR (%s)', behr_loc.ShortName), 'PECANS (\tau = 3h)');
            set(gca,'fontsize',16);
            
        end
        
        function [x, y] = center_and_normalize(x, y, varargin)
                % CENTER_AND_NORMALIZE Centers line densities on 0 and normalizes
                %
                %   [X, Y] = CENTER_AND_NORMALIZE(X, Y) Centers the line
                %   density define on coordinates X with values Y so that
                %   the max is at 0 and the values range from 0 to 1.
                %
                %   Parameters:
                %
                %       'norm' - changes normalization mode. Default is
                %       "squeeze", which compresses to the range 0 to 1.
                %       Other options:
                %           * 'max' - slides the line densities so that the
                %           max is at 0.
                %           * 'min' - slides the line densities so that the
                %           min is at 0.
                p = advInputParser;
                p.addParameter('norm', 'squeeze')
                p.parse(varargin{:});
                pout = p.Results;
                
                E = JLLErrors;
                
                norm_mode = pout.norm;
                allowed_norm_modes = {'squeeze', 'max', 'min'};
                if ~ismember(norm_mode, allowed_norm_modes)
                    E.badinput('"norm" must be one of "%s"', strjoin(allowed_norm_modes, '", "'));
                end
                
                [maxy,m] = max(y);
                x = x - x(m);
                switch lower(norm_mode)
                    case 'squeeze'
                        y = scale_to_range(y, [0 1]);
                    case 'max'
                        y = y - maxy;
                    case 'min'
                        y = y - min(y);
                    otherwise
                        E.notimplemented('No method defined for norm_mode = %s', norm_mode)
                end
            end
        
        function plot_each_pecans_fit()
            taus = 1:9;
            wind_speed_ms = 5;
            
            ld_fits_fig = figure;
            ld_fits_ax = gca;
            tau_v_tau_fig = figure;
            tau_v_tau_ax = gca;
            a_vs_tau_fig = figure;
            a_vs_tau_ax = gca;
            
            %pecans_results = repmat(make_empty_struct_from_cell({'x','linedens','ffit','emgfit','fit_stats'}), numel(taus), 1);
            pecans_fits_taus = nan(size(taus));
            pecans_fits_r2 = nan(size(taus));
            pecans_calc_a = nan(size(taus));
            line_colors = lines(); % remember all the color map strings correspond to functions that return the necessary matrix
            l = gobjects(numel(taus),1);
            legstr = cell(1, numel(taus));
            
            for i_tau = 1:numel(taus)
                pecans_data = misc_pecans_lifetime_plots.read_pecans_for_tau(taus(i_tau));
                [ffit, emgfit, fit_stats] = fit_line_density(pecans_data.x, pecans_data.A, 'none'); 
                
                % Cut down the fit to just where A is > 0.1% of max
                xx = pecans_data.A > 0.001 * max(pecans_data.A);
                [ffit_sub, emgfit_sub, fit_stats_sub] = fit_line_density(pecans_data.x(xx), pecans_data.A(xx), 'none'); 
                
                pecans_fits_taus(i_tau) = ffit.x_0 * 1000 / (wind_speed_ms * 3600);
                pecans_fits_r2(i_tau) = fit_stats.r2;
                pecans_calc_a(i_tau) = trapz(pecans_data.x, pecans_data.A);
                
                line(ld_fits_ax, pecans_data.x, pecans_data.A, 'color', line_colors(i_tau, :), 'linestyle', 'none', 'marker', 'o');
                l(i_tau) = line(ld_fits_ax, pecans_data.x, emgfit, 'color', line_colors(i_tau, :), 'linestyle', '-');
                line(ld_fits_ax, pecans_data.x(xx), emgfit_sub, 'color', line_colors(i_tau,:), 'linestyle', '--');
                legstr{i_tau} = sprintf('\tau = %d', taus(i_tau));
            end
            
            scatter(tau_v_tau_ax, taus, pecans_fits_taus, [], pecans_fits_r2,'filled');
            colormap(tau_v_tau_ax,'cool')
            cb = colorbar(tau_v_tau_ax);
            cb.Label.String = 'R2';
            xlabel('Model \tau');
            ylabel('Fit \tau');
            set(tau_v_tau_ax,'fontsize',16)
            title(tau_v_tau_ax, 'Fitted \tau vs. model \tau using all points')
            
            
            scatter(a_vs_tau_ax, taus, pecans_calc_a);
            plot_fit_line(a_vs_tau_ax, 'one2one', false);
            xlabel('Model \tau');
            ylabel('Calculated a');
            set(a_vs_tau_ax,'fontsize',16)
            title(a_vs_tau_ax, 'Calculated a vs model \tau')
        end
        
        function test_fitting_windows
            all_windows = [20, 30:2:50, 60:10:90];
            actual_n_points = nan(numel(all_windows),1);
            results = nan(numel(all_windows),4);
            keep = true(size(actual_n_points));
            for i_win = 1:numel(all_windows)
                window = all_windows(i_win);
                fprintf('Working on %d wide windows\n', window);
                wstate = warning('off');
                fits = misc_pecans_lifetime_plots.test_fit_quality_criteria('do_plots',false,'n_downwind_points',window,'allow_any_n_pts',true);
                warning(wstate);
                fits_accepted = [fits.is_fit_good];
                fits_correct = [fits.is_fit_correct];
                
                actual_n_points(i_win) = numel(fits(1).x);
                if i_win > 1 && actual_n_points(i_win) == actual_n_points(i_win-1)
                    fprintf('%d wide downwind window resulted in the same number of actual points used as the previous one, aborting further windows\n', window);
                    keep(i_win:end) = false;
                    break
                end
                % number of false positives, when the critera say the fit
                % is acceptable, but the lifetime and maxima are wrong.
                results(i_win,1) = sum(fits_accepted & ~fits_correct);
                % number of false negatives, when the critera say the fit
                % is not acceptable, but the lifetime and maxima are right
                results(i_win,2) = sum(~fits_accepted & fits_correct);
                % number of true negatives and true positives.
                results(i_win,3) = sum(~fits_accepted & ~fits_correct);
                results(i_win,4) = sum(fits_accepted & fits_correct);
            end
            
            figure;
            bar(actual_n_points(keep), results(keep,:), 'stacked');
            % scale the y axis so that the top limit is the number of
            % simulations
            xlabel('Window width');
            legend('False positives','False negatives','True negatives','True positives');
        end
        
        function fit_vars = test_fit_quality_criteria(varargin)
            E = JLLErrors;
            p = advInputParser;
            p.addParameter('do_plots', true);
            % want 60-70 total points, the cutting down algorithm adds 1/2
            % this number to upwind, so this will give us 61 total
            % (counting the center point)
            p.addParameter('n_downwind_points', 40);
            p.addParameter('allow_any_n_pts', false);
            p.parse(varargin{:});
            pout = p.Results;
            
            n_downwind_points = pout.n_downwind_points;
            do_plots = pout.do_plots;
            allow_any_n_pts = pout.allow_any_n_pts;
            
            % Load each of the lifetime model runs, record the line
            % density, x-coordinate, lifetime, and emission width of each.
            F = dirff(fullfile(misc_pecans_lifetime_plots.workspace_dir, 'PECANS', 'lifetime-ensemble', '*.nc'));
            fit_vars = repmat(make_empty_struct_from_cell({'x', 'A', 'tau', 'emission_width', 'fit_info', 'is_fit_good', 'is_lifetime_correct'}), numel(F), 1);
            for i_file = 1:numel(F)
                data = misc_pecans_lifetime_plots.read_pecans(F(i_file).name);
                
                % Cut down the model domain to the standard 100 km
                % upwind/200 km downwind box for cities.
                x = data.x;
                linedens = data.A;
                [~, i_max] = max(linedens);
                xx_box = floor(i_max - 0.5*n_downwind_points):ceil(i_max + n_downwind_points);
                xx_box = xx_box(xx_box >= 1 & xx_box <= length(x));
                %x_max = x(i_max);
                %xx_box = x >= x_max - 200e3 & x <= x_max + 400e3;
                
                % Add some noise to the line densities to mimic the fitting
                % that the algorithm actually has to do. Set the seed so
                % that the noise is reproducible each time.
                %
                % TODO: make sure this is a reasonable amount of noise that
                % does a good job mimicking the noise of the actual data
                rng(0);
                noise_frac = 0.05;
                % Allow positive or negative variation by noise_frac.
                linedens = linedens .* (1 + 2 * noise_frac * rand(size(linedens)) - noise_frac);
                
                fit_vars(i_file).x = x(xx_box) / 1000; % convert meters -> kilometers
                fit_vars(i_file).A = linedens(xx_box);
                fit_vars(i_file).tau = data.CONFIG.CHEMISTRY.mechanism_opts.lifetime_seconds / 3600; % convert seconds -> hours
                fit_vars(i_file).emission_width = data.CONFIG.EMISSIONS.emission_opts.width_x / 1000; % convert meters -> kilometers
                
                [fit_info.ffit, fit_info.emgfit, fit_info.param_stats, fit_info.f0, fit_info.history, fit_info.fitresults] = fit_line_density(fit_vars(i_file).x, fit_vars(i_file).A, 'none');
                fit_vars(i_file).fit_info = fit_info;
                fit_vars(i_file).is_fit_good = misc_emissions_analysis.is_fit_good(fit_vars(i_file).x, fit_vars(i_file).A, fit_vars(i_file).fit_info, 'any_num_pts', allow_any_n_pts);
                
                wind_speed = data.CONFIG.TRANSPORT.wind_speeds.x;
                fit_tau = fit_info.ffit.x_0 * 1000 / wind_speed / 3600; % wind speed in m/s, x_0 in, need to convert to h
                % Check that the lifetimes are within 10%
                fit_vars(i_file).is_lifetime_correct = abs(reldiff(fit_tau, fit_vars(i_file).tau)) < 0.1;
                % Check that the fit and line density maxima are within 2
                % points
                [~, ld_imax] = max(linedens(xx_box));
                [~, fit_imax] = max(fit_info.emgfit);
                fit_vars(i_file).is_max_correct = abs(ld_imax - fit_imax) <= 2;
                
                fit_vars(i_file).is_fit_correct = fit_vars(i_file).is_lifetime_correct && fit_vars(i_file).is_max_correct;
                fit_vars(i_file).fit_info.tau = fit_tau;
            end
            
            % Okay, so now that we've read all the files, we need to list
            % all of the lifetimes and emissions widths so that we can make
            % a grid plot showing cases where the fit fails and is
            % rejected, where the fit fails and is not rejected, where the
            % fit works and is correctly kept, and where the fit works and
            % is incorrectly rejected.
            
            if do_plots
                plot_data = struct('correctly_kept', struct('color',[0 0.5 0], 'taus', [], 'emwidths', [], 'fit_taus', []),...
                    'incorrectly_rejected', struct('color', [1 0.5 0], 'taus', [], 'emwidths', [], 'fit_taus', []),...
                    'correctly_rejected', struct('color', 'c', 'taus', [], 'emwidths', [], 'fit_taus', []),...
                    'incorrectly_kept', struct('color', 'y', 'taus', [], 'emwidths', [], 'fit_taus', []));
                
                for i_sim = 1:numel(fit_vars)
                    this_fit_good = fit_vars(i_sim).is_fit_good;
                    this_fit_correct = fit_vars(i_sim).is_fit_correct;
                    
                    if this_fit_good && this_fit_correct
                        field = 'correctly_kept';
                    elseif ~this_fit_good && this_fit_correct
                        field = 'incorrectly_rejected';
                    elseif ~this_fit_good && ~this_fit_correct
                        field = 'correctly_rejected';
                    elseif this_fit_good && ~this_fit_correct
                        field = 'incorrectly_kept';
                    else
                        E.notimplemented('Should not have reached here')
                    end
                    
                    plot_data.(field).taus(end+1) = fit_vars(i_sim).tau;
                    plot_data.(field).emwidths(end+1) = fit_vars(i_sim).emission_width;
                    plot_data.(field).fit_taus(end+1) = fit_vars(i_sim).fit_info.tau;
                end
                
                figure;
                grid_ax = gca;
                fns = fieldnames(plot_data);
                l = gobjects(numel(fns),1);
                legstr = cell(1, numel(fns));
                keep = true(size(l));
                
                tau_fig = figure;
                for i_field = 1:numel(fns)
                    figure(tau_fig);
                    tau_ax = subplot(2,2,i_field);
                    this_field = fns{i_field};
                    if isempty(plot_data.(this_field).taus)
                        keep(i_field) = false;
                        continue
                    end
                    l(i_field) = line(grid_ax, plot_data.(this_field).taus, plot_data.(this_field).emwidths, 'color', plot_data.(this_field).color, 'markerfacecolor', plot_data.(this_field).color, 'marker', 'o','linestyle','none');
                    legstr{i_field} = capitalize_words(strrep(this_field, '_', ' '));
                    
                    scatter(tau_ax, plot_data.(this_field).taus, plot_data.(this_field).fit_taus);
                    xlabel(tau_ax, 'Model \tau (h)');
                    ylabel(tau_ax, 'Fit \tau (h)');
                    title(legstr{i_field});
                end
                
                legend(grid_ax, l(keep), legstr(keep));
                xlabel(grid_ax, 'Model \tau (h)')
                ylabel(grid_ax, 'Model emission width (\sigma, km)')
            end
        end
        
        function fit_test = plot_fit_test(fit_test)
            if ~isscalar(fit_test)
                taus = [fit_test.tau];
                emwidths = [fit_test.emission_width];
                u_taus = cellfun(@num2str, num2cell(unique(taus)), 'uniform', false);
                u_emwidths = cellfun(@num2str, num2cell(unique(emwidths)), 'uniform', false);
                
                if numel(u_taus) > 1
                    desired_tau = str2double(ask_multichoice('Which model lifetime?', u_taus, 'list', true));
                else
                    desired_tau = str2double(u_taus{1});
                end
                
                if numel(u_emwidths) > 1
                    desired_emwidth = str2double(ask_multichoice('Which emission width?', u_emwidths, 'list', true));
                else
                    desired_emwidth = str2double(u_emwidths{1});
                end
                
                xx = taus == desired_tau & emwidths == desired_emwidth;
                fit_test = fit_test(xx);
            end
            figure; 
            plot(fit_test.x, fit_test.A, 'bo', fit_test.x, fit_test.fit_info.emgfit, 'r-')
        end
        
        %%%%%%%%%%%%%%%%%%%%
        % Helper functions %
        %%%%%%%%%%%%%%%%%%%%
        
        function data = read_pecans_for_tau(tau, varargin)
            p = advInputParser;
            p.addOptional('emwidth', 3.0);
            p.parse(varargin{:});
            pout = p.Results;
            emwidth = pout.emwidth;
            pecans_name = sprintf('pecans_ens_tau-%.1fh_emwidth-%.1fkm.nc', tau, emwidth);
            data = misc_pecans_lifetime_plots.read_pecans(fullfile(misc_pecans_lifetime_plots.pecans_dir, 'lifetime-ensemble', pecans_name));
            % convert x from meters to kilometers
            data.x = data.x/1000;
        end
        
        function data = read_pecans(pecans_file)
            pinfo = ncinfo(pecans_file);
            for i_var = 1:numel(pinfo.Variables)
                var_name = pinfo.Variables(i_var).Name;
                data.(var_name) = ncread(pinfo.Filename, var_name);
            end
            data.CONFIG = misc_pecans_lifetime_plots.parse_pecans_config(pecans_file);
        end
        
        function PecansConfig = parse_pecans_config(pecans_file)
            config_string = strsplit(ncreadatt(pecans_file, '/', 'pecans_config'), '\n');
            PecansConfig = struct();
            section_name = '';
            for i_line = 1:numel(config_string)
                this_line = config_string{i_line};
                % Is this a section header: "[SECTION]"
                if regcmp(this_line, '\[.+\]')
                    if ~isempty(section_name)
                        PecansConfig.(section_name) = current_section;
                    end
                    
                    section_name = regexp(this_line, '(?<=\[).*(?=\])', 'match', 'once');
                    current_section = struct();
                % If not, then it must be an option line (if not empty)
                elseif ~isempty(this_line)
                    split_line = strsplit(this_line, '=');
                    option_name = strtrim(split_line{1});
                    value = strtrim(split_line{2});
                    
                    % The value may be a string, a number, or a dictionary.
                    % Dictionaries are enclosed in {}. 
                    if regcmp(value, '\{.+\}')
                        current_section.(option_name) = parse_dict(value);
                    else
                        current_section.(option_name) = parse_value(value);
                    end
                end
            end
            
            % Make sure we include the last section
            PecansConfig.(section_name) = current_section;
            
            
            function val = parse_value(str_in)
                % Try to convert the string to a number. If it cannot be
                % parsed as a number, it must be a string.
                val = str2double(str_in);
                if isnan(val) 
                    val = str_in;
                end
            end
            
            function dict = parse_dict(str_in)
                dict = struct();
                str_in = regexp(str_in, '(?<=\{).*(?=\})', 'match', 'once');
                keyval_pairs = strsplit(str_in, ',');
                for i_pair = 1:numel(keyval_pairs)
                    split_pair = strsplit(keyval_pairs{i_pair}, ':');
                    % Remove any extra quotes, either single or double
                    key = regexp(strtrim(split_pair{1}), '(?<=[''"]).*(?=[''"])', 'match', 'once');
                    val = strtrim(split_pair{2});
                    dict.(key) = parse_value(val);
                end
            end
        end
    end
end

