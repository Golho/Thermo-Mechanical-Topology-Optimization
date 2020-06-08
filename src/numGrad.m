function [grad] = numGrad(func, x, h)
%UNTITLED Summary of this function goes here
%   Detailed explanation goes here
grad = zeros(size(x));
for i = 1:numel(x)
    dx = x;
    dx(i) = dx(i) + h/2;
    f1 = func(dx);
    dx(i) = dx(i) - h;
    f2 = func(dx);
    grad(i) = (f1 - f2) / h;
end
end

