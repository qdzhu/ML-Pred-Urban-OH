function [outputArg1,outputArg2] = emission_analysis_plots(varargin)
%UNTITLED Summary of this function goes here
%   Detailed explanation goes here

p = advInputParser;
p.addFlag('remake_data');
p.addFlag('overwrite_data');

p.parse(varargin{:});
pout = p.AdvResults;

remake_data = pout.remake_data;
do_overwrite = pout.overwrite_data;


if remake_data
    run_local_data_prep(do_overwrite);
end

    function run_local_data_prep(do_overwrite)
        
        wrf_bools = {true, false};
        loc_inds = {9, 1:70}; % 26 Feb 2018, I only have WRF line densities for Chicago
        
        for i_wrf = 1:numel(wrf_bools)
            time_opts = {'beginning','UMTWRFS';...
                'beginning', 'TWRF';...
                'beginning', 'US';...
                'end', 'UMTWRFS';...
                'end', 'TWRF';...
                'end', 'US'};
            n_times = size(time_opts,1);
            
            options = struct('time_period', time_opts(:,1),...
                'use_wrf', repmat(wrf_bools(i_wrf),n_times,1),...
                'loc_indicies', repmat(loc_inds(i_wrf),n_times,1),...
                'add_nei', repmat({true},n_times,1),...
                'days_of_week', time_opts(:,2),...
                'do_overwrite', repmat({do_overwrite}, n_times, 1),...
                'fatal_fit_fail', repmat({false}, n_times, 1));
            
            for i_opts = 1:numel(options)
                try
                    misc_emissions_analysis.make_emg_fits(options(i_opts));
                catch err
                    if strcmp(err.identifier, 'emis_analysis:no_linedens_file')
                        fprintf('%s\n', err.message);
                        continue
                    else
                        rethrow(err);
                    end
                end
            end
        end
    end



end

