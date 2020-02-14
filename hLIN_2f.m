function [fe] = hLIN_2f(ex, ey, loadMagnitude, t)
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
    

I = 1/2 * [1;
           1];
fe = t * length * loadMagnitude * I;
end
