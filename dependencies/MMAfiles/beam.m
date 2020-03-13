%    This is the file beam.m
%    which defines the cantilever beam problem.
%
function [f0val,df0dx,df0dx2,fval,dfdx,dfdx2] = beam(x);
%
%   written in January 1999 by
%
%   Krister Svanberg (krille@math.kth.se)
%   Optimization and Systems Theory, KTH,
%   SE-10044 Stockholm, Sweden.
%
e = [ 1 1 1 1 1 ]';
df0dx = e;
df0dx2 = 0*e;
f0val = df0dx'*x;
coef = [61 37 19 7 1]';
x2 = x.*x;
x3 = x2.*x;
x4 = x2.*x2;
x5 = x3.*x2;
x3inv = e./x3;
fval = coef'*x3inv - 1;
dfdx = -3*(coef./x4)';
dfdx2 = 12*(coef./x5)';


