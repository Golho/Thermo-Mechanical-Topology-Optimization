function [Ke] = flw2re(ex, ey, D, t)
%FLW2RE Compute the element heat stiffness matrix for a 2D Melosh element
% It is assumed that the conductivity matrix is symmetric.
%   ex      x-coordinates of the nodes [m]
%   ey      y-coordinates of the nodes [m]
%   D       Thermal conductivity matrix [2 x 2]
%   t       Thickness [m]
area = abs(ex(3)-ex(1))*abs(ey(3)-ey(1));

s1 = D(2, 2)/6  -    D(1, 1)/3;
s2 = -D(2, 2)/3 +    D(1, 1)/6;
s3 = -D(2, 2)/6 -    D(1, 1)/6      -   D(1, 2)/2;
s4 = -D(2, 2)/6 -    D(1, 1)/6      +   D(1, 2)/2;
s5 = D(2, 2)/3  +    D(1, 1)/3      -   D(1, 2)/2;
s6 = D(2, 2)/3  +    D(1, 1)/3      +   D(1, 2)/2;

I = [s6,     s1,     s3,     s2;
     s1,     s5      s2      s4;
     s3      s2      s6      s1;
     s2      s4      s1      s5];
Ke = t*area*I;
end
