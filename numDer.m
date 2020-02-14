function [der] = numDer(func, x)
%NUMDER Summary of this function goes here
%   Detailed explanation goes here
h = 0.01;
f1 = func(x-h/2);
f2 = func(x+h/2);
der = (f2-f1)/h;
end

