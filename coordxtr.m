function [Ex, Ey, Ez] = coordxtr(Enod, Coord)
%COORDXTR Summary of this function goes here
%   Detailed explanation goes here
spatialDimension = size(Coord, 2);
xCoord = Coord(:, 1);
Ex = xCoord(Enod(:, 2:end));
if spatialDimension >= 2
    yCoord = Coord(:, 2);
    Ey = yCoord(Enod(:, 2:end));
end
if spatialDimension == 3
    zCoord = Coord(:, 1);
    Ez = zCoord(Enod(:, 2:end));
end
end

