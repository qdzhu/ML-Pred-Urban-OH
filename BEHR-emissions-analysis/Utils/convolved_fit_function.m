function emgfxn = convolved_fit_function(x_slow, slow_ld)
%CONVOLVED_FIT_FUNCTION Create a function handled represented a convolved line density fit
%   EMGFIT = CONVOLVED_FIT_FUNCTION( X_SLOW, SLOW_LD ) given an NO2 line
%   density, SLOW_LD defined at x-coordinates X_SLOW, this function
%   generates a function handle that can be used as the fitting function
%   for a line density taken at fast winds. The returned function handle,
%   EMGFIT, takes two arguments F and X. F is the vector of five fitting
%   parameters (a, x_0, mu_x, sigma_x, and B) varied by FIT_LINE_DENSITY().
%   (sigma_x is retained for compatibility even though this function does
%   not actually use sigma_x.) X is the x-coordinates that the fast line
%   density being fit is defined on; it must be the same as the coordinates
%   that the slow line density was defined on or an error is thrown.

emgfxn = @convolved_fxn;

    function emg = convolved_fxn(f,x)
        if ~isequal(x, x_slow)
            error('convolved_emg_fit:x_mismatch', 'New x coordinates do not match those passed for the slow line densities')
        end
        
        a = f(1);
        x_0 = f(2);
        mu_x = f(3);
        B = f(5);
        
        exp_component = zeros(size(x));
        exp_component(x > mu_x) = exp( -(x(x>mu_x)-mu_x)/x_0 );
        emg = a ./ x_0 .* conv_trapz(x, exp_component, slow_ld) + B;
    end

end

