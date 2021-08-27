function [ emis, uncert_emis, tau, uncert_tau ] = compute_emg_emis_tau( a, uncert_a, x0, uncert_x0, wind_mode, varargin )
%[ EMIS, UNCERT_EMIS, TAU, UNCERT_TAU ] = COMPUTE_EMG_EMIS_TAU( A, UNCERT_A, X0, UNCERT_X0, 'vec', WIND_SPEED_VEC )
%   This function computes the values of emissions, lifetime, and their
%   uncertainties from the values of a, x_0, their uncertainties, and the
%   vector of wind speeds for days used in generating the line densities
%   that fit. For example, if windvel is a 90x1 vector of wind speeds (in
%   m/s) for the 90 days considered for the line densities and only days
%   with wind 3 m/s were used, then pass windvel(windvel>3) as
%   WIND_SPEED_VEC. This is used to compute the average wind speed and
%   error in the wind (as a 95% confidence interval) for use in the
%   uncertainty calculations. EMIS and UNCERT_EMIS are given in Mg/hr and
%   TAU and UNCERT_TAU are given in hours. WIND_SPEED_VEC must be given in
%   meters/second.
%
%[ EMIS, UNCERT_EMIS, TAU, UNCERT_TAU ] = COMPUTE_EMG_EMIS_TAU( A, UNCERT_A, X0, UNCERT_X0, 'avg', WIND_SPEED_MEAN, WIND_SPEED_ERROR)
%   In this format, the mean and error of the wind speed is given instead.
%   This can be useful if you just stored these values rather than the full
%   vector, or if you wish to use an alternate definition of error in the
%   wind speed. Both WIND_SPEED_MEAN and WIND_SPEED_ERROR must be in
%   meters/second.

%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%% INPUT CHECKING %%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%
E = JLLErrors;

p = inputParser;
p.addOptional('wind_vec_or_mean',[],@isnumeric);
p.addOptional('wind_error',[],@isnumeric);
p.addParameter('emissions_type','nox');

p.parse(varargin{:});
pout = p.Results;

wind_vec_or_mean = pout.wind_vec_or_mean;
wind_error = pout.wind_error;
emissions_type = pout.emissions_type;

if ~isequal(size(a), size(uncert_a)) || ~isequal(size(a),size(x0)) || ~isequal(size(a),size(uncert_x0))
    E.badinput('A, UNCERT_A, X0, UNCERT_X0 must all be the same size.')
end

allowed_wind_modes = {'vec','avg'};
if ~ismember(wind_mode, allowed_wind_modes)
    E.badinput('WIND_MODE (5th input) must be one of %s', strjoin(allowed_wind_modes, ', '));
end

switch lower(wind_mode)
    case 'vec'
        if isempty(wind_vec_or_mean)
            E.badinput('When using WIND_MODE == ''vec'' there must be one additional input, the vector of wind speeds');
        end
    case 'avg'
        if isempty(wind_vec_or_mean) || isempty(wind_error)
            E.badinput('When using WIND_MODE == ''avg'' there must be two additional inputs, the mean and error of wind speeds');
        elseif any(~iscellcontents({wind_vec_or_mean, wind_error}, 'isscalar'))
            E.badinput('When using WIND_MODE == ''avg'' both additional inputs are expected to be scalars.')
        end
end

%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%% MAIN FUNCTION %%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%

if strcmpi(wind_mode, 'vec')
    wind_speed_vec = wind_vec_or_mean;
    avg_wind = nanmean(wind_speed_vec/1000*3600); % windvel in m/s, convert to km/h
    student_t = tinv(0.975, numel(wind_speed_vec)); %tinv gives one-tailed, we want two-tailed 95% CI
    err_wind = (student_t * nanstd(wind_speed_vec/1000*3600))/sqrt(numel(wind_speed_vec));
else
    avg_wind = wind_vec_or_mean/1000*3600;
    err_wind = wind_error/1000*3600;
end

% Always convert NO2 to NOx - we want total emissions, whether it is given
% in moles, Mg NO-equivalent, Mg NO2-equivalent, or Mg of NOx.
emis = 1.32 .* a .* avg_wind ./ x0;  

if strcmpi(emissions_type, 'nox')
    % Calculate assumed mass of NOx if NOx:NO2 ratio is 1.32
    % MM NO = 30.01 g/mol
    % MM NO2 = 46.01 g/mol
    mol2Mg = (1/1.32 * 46.01 + (1-1/1.32)*30.01)*1e-6;
elseif strcmpi(emissions_type, 'no')
    mol2Mg = 30.01e-6;
elseif strcmpi(emissions_type, 'no2')
    mol2Mg = 46.01e-6;
end
emis = emis * mol2Mg;

tau = x0 ./ avg_wind;
% Uncertainty in emissions needs to add the uncertainty in the
% NOx:NO2 ratio (10%). Uncertainty in lifetime will just depend
% on uncertainty in x_0

%Simple case (assumes average wind has no error)
%uncert_tau = uncert_x0  ./ avg_wind;
%Full case (consider uncertainty as 95% CI
uncert_tau = sqrt( (uncert_x0 ./ avg_wind).^2 + ( err_wind .* -x0 ./ avg_wind.^2 ).^2 );

% We'll handle the uncertainty in E in two steps. First,
% compute the percent error due to error in a and tau. Then add
% this in quadrature with the 10% error due to the NOx:NO2
% ratio.
% init_uncert_E_squared = (uncert_a .* mol2Mg ./ tau).^2 + (uncert_tau .* -a .* mol2Mg ./ (tau .^ 2)).^2;
% per_uncert_E = sqrt( init_uncert_E_squared ./ (E .^ 2) + 0.1 .^ 2 );
% uncert_E = E .* per_uncert_E;
NOxNO2 = 1.32;
uncert_NOxNO2 = 0.1 * NOxNO2;
uncert_emis = sqrt((a .* mol2Mg ./ tau .* uncert_NOxNO2).^2 + (NOxNO2 ./ tau .* uncert_a .* mol2Mg).^2 + (-NOxNO2 .* a .* mol2Mg ./ tau.^2 .* uncert_tau).^2);

end

