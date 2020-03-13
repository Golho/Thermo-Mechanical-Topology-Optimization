function [Ke] = hTRI_3K(ex, ey, D, t)
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
if area <= 0
    error('The element has no area or has inverse node order');
end
    
% B-matrix for a triangular element with 3 nodes
Be = [ ey(2)-ey(3)     ey(3)-ey(1)      ey(1)-ey(2); 
       ex(3)-ex(2)     ex(1)-ex(3)      ex(2)-ex(1) ] * 1/(2*area);

I = Be' * D * Be;
Ke = t*area*I;
end
