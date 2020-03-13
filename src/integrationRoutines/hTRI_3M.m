function [Me] = hTRI_3M(ex, ey, rho, t)
%FLW2RE Compute the element heat stiffness matrix for a 2D Melosh element
% It is assumed that the conductivity matrix is symmetric.
%   ex      x-coordinates of the nodes [m]
%   ey      y-coordinates of the nodes [m]
%   D       Thermal conductivity matrix [2 x 2]
%   t       Thickness [m]

% reassure the correct dimensions
ex = reshape(ex, 3, 1);
ey = reshape(ey, 3, 1);

% area of the element
area = 1/2*det([ones(3,1) ex ey]);
    

I = 1/12 * [2 1 1;
            1 2 1;
            1 1 2];
Me = t * area * rho * I;
end
