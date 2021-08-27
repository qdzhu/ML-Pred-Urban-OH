function c = conv_trapz(x, f, g)
%CONV_TRAPZ Approximate convolution of two continuous functions using trapezoidal integration.
%   C = CONV_TRAPZ( X, F, G ) Convolves two functions, F and G, represented
%   as values on a discrete grid defined at coordinates X. Approximates the
%   convolution integral using the TRAPZ function and returns convolution
%   result as an array C the same size as X, F, and G. Unlike CONV(), this
%   does not have any options to include the zero-padded part of the
%   convolution. Note that internally X, F, and G are all converted to
%   column vectors,
%
%   There is a subtle but important difference between this function and
%   CONV() that is best explained by example. Consider the convolution of a
%   boxcar function with itself:
%
%       x = -2:0.5:2; 
%       f = double(abs(x) <= 0.5)
%
%   When thinking about convolution of two continuous functions, we would
%   expect the result to be a piecewise function that increases linearly
%   from 0 to a maximum then decreases back to 0, since the value of the
%   convolution at some position t represents the overlap between F and G
%   if G is reversed and shifted by t. So the maximum value comes when F
%   and G maximally overlap. If F and G are the same, and are a boxcar
%   function with width and height 1, then that maximum value should be 1
%   and it should be 0 at x < -1 and x > 1 (since when G is shifted by at
%   least the width of the boxcar there is no overlap). 
%
%   However, because CONV() does discrete convolution, then the maximum is
%   still when G isn't shifted at all (in this case), but its max value is
%   3, not 1, because it is calculated as sum(f .* f) and there are three
%   values of 1 and the rest are 0. Worse, if you gave f finer resolution,
%   e.g.
%
%       x = -2:0.01:2;
%       f = double(abs(f)<=0.5);
%
%   then CONV(f,f) would give a vector with a maximum of 101, simply
%   because there are more points. This may be useful if your function
%   truly is a series of discrete impulses, but gives counterintuitive
%   results if trying to apply it to a physically continuous function that
%   has been sampled at discrete points.
%
%   If you use the second x and f above and compare:
%
%       c1 = conv(f,f);
%       c2 = conv_trapz(x,f,f);
%
%   you'll see that both have the expected shape, but c1's maximum is much
%   too large, whereas conv_trapz is much closer to the correct maximum.
%   (It will be a little off because f is not a true boxcar function, it
%   doesn't jump at x = +/-0.5, there's a short ramp between -0.51 and
%   -0.5/0.5 and 0.51. Making f coarser makes the agreement worse.)

E = JLLErrors;
if ~isequal(size(x), size(f)) || ~isequal(size(x), size(g))
    E.badinput('X, F, and G must all be the same size')
elseif ~isvector(x) || ~isvector(f) || ~isvector(g)
    E.badinput('X, F, and G must all be vectors')
elseif ~isnumeric(x) || ~isnumeric(f) || ~isnumeric(g)
    E.badinput('X, F, and G must all be numeric')
end



% The definition of a convolution is:
%
%   (f * g)(t) = \int_{-\infty}^{+\infty} f(\tau) g(t - \tau) d\tau
%
% in other words, the valus of the convolution f * g at t is equal to the
% integral of f times g shifted by t and reversed.
%
% We handle the shifting of g by interpolating it to the shifted
% coordinates, but element-wise multiplying it with f as if it were still
% defined on the same x.

c = nan(size(x));
for i=1:numel(x)
    xg = x(i) - x;
    ginterp = interp1(x,g,xg,'linear',0);
    c(i) = trapz(x, f .* ginterp);
end

end

