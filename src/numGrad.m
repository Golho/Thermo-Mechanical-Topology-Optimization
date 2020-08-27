function [grad] = numGrad(func, x, h)
%NUMGRAD    Numerical gradient of a function
%   grad = numGrad(func, x, h) Numerically calculates the gradient of a 
%   multidimensional function, func, at a point x. The parameter h decides
%   the step taken to calculate the slope.

grad = zeros(size(x));
for i = 1:numel(x)
    dx = x;
    dx(i) = dx(i) + h/2;
    f1 = func(x + h/2);
    dx(i) = dx(i) - h;
    f2 = func(dx);
    grad(i) = (f1 - f2) / h;
end
end

