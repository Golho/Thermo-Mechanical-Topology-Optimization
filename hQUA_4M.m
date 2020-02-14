function [Me] = hQUA_4M(ex, ey, rho, t)
%FLW2RE Compute the element heat mass matrix M for a 2D Melosh element
%   ex      x-coordinates of the nodes [m]
%   ey      y-coordinates of the nodes [m]
%   rho     Density or equivalent in integral M = \int_V N^T rho N dV
%   t       Thickness [m]
area = abs(ex(3)-ex(1))*abs(ey(3)-ey(1));

I = [4 2 1 2;
     2 4 2 1;
     1 2 4 2;
     2 1 2 4] * 1/36;
Me = t * rho * area * I;
end
