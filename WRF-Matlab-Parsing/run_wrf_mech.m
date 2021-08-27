function [ conc_out, species  ] = run_wrf_mech( model_name, init_cond )
%UNTITLED10 Summary of this function goes here
%   Detailed explanation goes here

% Initial conditions. Concentrations should be in molec. cm^{-3}. Each row
% of the cell array should have the name of the species (case insensitive)
% and the concentration.

Nair = 2e19; % number density of air at surface
T = 298;

if ~exist('model_name','var')
    model_name = 'noxbox';
end
if ~exist('init_cond','var')
    init_cond = {   'O3', 40e-9 * Nair;...
                    'NO2', 0;...
                    'NO', 1e-9 * Nair;...
                    'M', Nair};
end

solver = 'ode15s'; % can be rk4, simple, or ode15s

[ J, species, isfixed, photo_calls ] = parse_wrf_mech(model_name); %#ok<ASGLU>

lon = -84.39;
lat = 33.775;
date_in = '2013-06-01';

nt = 1800; % number timesteps
dt = 1; % seconds
save_freq = 1; % how frequently to save out the concentrations
start_utc_time = 19; % hr utc
run_mode = 'normal';
%%%%%%%%%%%%%%%%%%
%%% Initialize %%%
%%%%%%%%%%%%%%%%%%

C = zeros(size(species)); % concentration vector
for a=1:size(init_cond,1)
    xx = strcmpi(init_cond{a,1},species);
    C(xx) = init_cond{a,2};
end

conc_out = nan(numel(C), ceil(nt/save_freq));

tuv_hr = nan;

%%%%%%%%%%%
%%% RUN %%%
%%%%%%%%%%%

% Make the photolysis a nested variable for persistence
j = 0;
switch run_mode
    case 'normal'
        % Run in regular mode
        for t=1:nt
            utc_time = start_utc_time + t*dt/(3600);
            mech_timestep;
        end
        
    case 'steady-state'
        % Run in SS mode
        E.notimplemented('%s','End condition for steady state model not implemented, would run forever');
        utc_time = start_utc_time;
        while true
            mech_timestep;
        end
end


function mech_timestep
    % Should we call TUV again? Do so if we are closer to a new hour than
    % we were before.
    if round(utc_time) ~= tuv_hr
        tuv_hr = round(utc_time);
        j = call_tuv(photo_calls, date_in, tuv_hr, lon, lat, true);
    end
    dC = zeros(size(C));
    for i=1:numel(C)
        switch lower(solver)
            case 'simple'
                dC(i) = J{i}(C,j,T,Nair)*dt;
            case 'rk4'
                dC_dt = @(t_n, c_n) J{i}(c_n, j, T, Nair);
                dC(i) = int_rk4(dC_dt, dt, 0, C); % change in concentration does not directly depend on time, so pass a dummy value
            case 'ode15s'
                dC_dt = @(t, c) ode15s_fun(J, c, j, T, Nair);
                [tout, cout] = ode15s(dC_dt, [0 dt], C);
                dC = cout(end,:) - cout(1,:); % necessary to be consistent with other methods.
                break % only need to call this once to get all species.
        end
    end
    C = C + dC;
    
    % Save if requested
    if mod(t,save_freq) == 0
        conc_out(:,ceil(t/save_freq)) = C;
    end
end

end

function dcdt = ode15s_fun(J, c, j, T, Nair)
% Sub function used to generate the derivatives in an intermediate format
% that can then be wrapped in an anonymous function to pass to ode15s as
% dcdt = @(t, c) ode15s_fun(J, c, j, T, Nair)
dcdt = zeros(size(J));
colbool = false;
if iscolumn(c); c = c'; colbool = true; end
for i=1:numel(J)
    dcdt(i) = J{i}(c, j, T, Nair);
end
if colbool && isrow(dcdt); dcdt = dcdt'; end 

end