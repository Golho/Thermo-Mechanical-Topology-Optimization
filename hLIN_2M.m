function [Me] = hLIN_2M(ex, ey, rho, t)
%FLW2RE Compute the element heat stiffness matrix for a 2D Melosh element
% It is assumed that the conductivity matrix is symmetric.
%   ex      x-coordinates of the nodes [m]
%   ey      y-coordinates of the nodes [m]
%   rho     Density or equivalent
%   t       Thickness [m]

% reassure the correct dimensions
ex = reshape(ex, 2, 1);
ey = reshape(ey, 2, 1);

% length of the element
length = sqrt(diff(ex)^2 + diff(ey)^2);
    

I = 1/6 * [2 1;
           1 2];
Me = t * length * rho * I;
end
