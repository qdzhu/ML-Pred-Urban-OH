function [ fxn ] = wrf_rate_expr( law_name, varargin )
%FXN = WRF_RATE_EXPR( LAW_NAME )  Return a function handle for WRF-Chem KPP rate constants.
%   WRF-Chem's KPP defines non-photolytic rate constants using one of a
%   number of built-in functions to compute the numeric rate constant based
%   on physical parameters defined in the .eqn file. This form of this
%   function will return a function handle to the requested rate function.
%   LAW_NAME is one of the strings: 'arr2', 'k46', 'ko1d', 'troe', or
%   'troee' returning the correspondingly named function.
%
%   K = WRF_RATE_EXPR( LAW_NAME, PARAMETERS ) will evaluate the rate
%   constant based on the values of the PARAMETERS specified. Each rate law
%   has a specific set of parameters it required. You can call the help
%   text for any of the rate law functions with help
%   wrf_rate_expr>law_name.

switch lower(law_name)
    case 'arr2'
        fxn = @ARR2;
    case 'k46'
        fxn = @k46;
    case 'ko1d'
        fxn = @ko1d;
    case 'troe'
        fxn = @TROE;
    case 'troee'
        fxn = @TROEE;
    case 'fall'
        fxn = @FALL;
    otherwise
        error('wrf_rate_expr:badinput','%s not a valid rate expression name',law_name)
end

if numel(varargin) > 0
    fxn = fxn(varargin{:});
end

end

% These are defined in the .def file.

function k = k46(TEMP, C_M)
% k = k46(temperature, cair)
%   A fairly complex rate constant. Temperature must be in K, cair is the
%   number density of air in molec. cm^{-3}.
if nargin < 2
    error('rate_law:not_enough_inputs','k46 requires two inputs: TEMP, C_M (number density of air in molec./cm^3)')
end
k0=2.4e-14 .* exp(460 ./ TEMP);
k2=2.7e-17 .* exp(2199 ./ TEMP);
k3=6.5e-34 .* exp(1335 ./ TEMP) .* C_M;

k=k0+k3/(1+k3/k2);
end

function k = ko1d(TEMP, C_M)
% k = ko1d(temperature, cair)
%   Used solely for the relaxation of O(1D) to O(3P). Temperature must be
%   in K, cair is the number density of air in molec. cm^{-3}.
if nargin < 2
    error('rate_law:not_enough_inputs','ko1d requires two inputs: TEMP, C_M (number density of air in molec./cm^3)')
end
kN = 0.78084 .* C_M .* 1.8e-11 * exp(107 ./ TEMP);
k0 = 0.20946 .* C_M .* 3.2e-11 * exp(67 ./ TEMP);
k = kN + k0;
end

% These are defined in the WRFUserRateLaws.f90 file (in
% WRFV3/chem/KPP/kpp/kpp-2.1/util/WRF_Conform). 
function k = ARR2(A0, B0, TEMP)
% k = ARR2(A0, B0, temperature)
%    A simplified Arrhenius expression; k = A0 * exp(B0 / TEMP). 
if nargin < 3
    error('rate_law:not_enough_inputs','ARR2 requires three inputs: A0, B0, TEMP')
end
k = A0 * exp(-B0 ./ TEMP);
end

function k = TROE(k0_300K, n, kinf_300K, m, temp, cair)
% k = TROE(k0_300K, n, kinf_300K, m, temperature, cair)
%   Very close to the JPL definition for a 3-body reaction, i.e. k0_T =
%   k0_300K * (T/300)^(-n) and kinf_T = kinf_300K * (T/300)^(-m),
%   representing the low and high pressure limits. These are combined to
%   give the final rate constant. Temperature must be in K, cair is the
%   number density of air in molec. cm^{-3}.
if nargin < 6
    error('rate_law:not_enough_inputs','TROE requires 6 inputs: k0_300K, n, kinf_300K, m, temp, cair (number density of air in molec./cm^3)')
end
zt_help = 300 ./ temp;
k0_T    = k0_300K   .* zt_help .^ (n) .* cair; % k_0   at current T
kinf_T  = kinf_300K .* zt_help .^ (m);     % k_inf at current T
k_ratio = k0_T ./ kinf_T;
k   = k0_T ./ (1 + k_ratio) .* 0.6 .^ (1 ./ (1+log10(k_ratio).^2));
end

function k = TROEE(A, B, k0_300K,n,kinf_300K,m,temp,cair)
% k = TROE(A, B, k0_300K, n, kinf_300K, m, temperature, cair)
%   Similar to TROE, but for an equilibrium. k0_T = k0_300K * (T/300)^(-n)
%   and kinf_T = kinf_300K * (T/300)^(-m), representing the low and high
%   pressure limits. These are combined to give the final rate constant.
%   Temperature must be in K, cair is the number density of air in molec.
%   cm^{-3}.
if nargin < 8
    error('rate_law:not_enough_inputs','TROEE requires 8 inputs: A, B, k0_300K, n, kinf_300K, m, temp, cair (number density of air in molec./cm^3)')
end
zt_help = 300 ./ temp;
k0_T    = k0_300K   .* zt_help .^ (n) .* cair; % k_0   at current T
kinf_T  = kinf_300K .* zt_help .^ (m);       % k_inf at current T
k_ratio = k0_T ./ kinf_T;
troe   = k0_T ./ (1+k_ratio).*0.6.^(1 ./ (1+log10(k_ratio).^2));

k = A .* exp( - B ./ temp) .* troe;
end

function k = FALL(A0, B0, C0, A1, B1, C1, CF, temp, cair)
if nargin < 9
    error('rate_law:not_enough_inputs','FALL requires 8 inputs: A0, B0, C0, A1, B1, C1, CF, temp, cair (number density of air in molec./cm^3)')
end
    k0 = A0 .* exp(-B0/temp) .* (temp ./ 300).^C0;
    k1 = A1 .* exp(-B1/temp) .* (temp ./ 300).^C1;
    k0 = k0 .* cair;
    kratio = k0 ./ k1;
    k = (k0 ./ (1 + kratio)) .* CF .^ (1 ./ (1 + log10(kratio).^2));
end
