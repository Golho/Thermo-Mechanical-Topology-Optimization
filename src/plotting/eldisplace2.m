function [plot] = eldisplace2(Ex, Ey, ed, scaleFactor)
%ELDISPLACEMENT2 Summary of this function goes here
%   Detailed explanation goes here
Ex = Ex + scaleFactor*ed(1:2:end, :);
Ey = Ey + scaleFactor*ed(2:2:end, :);

plot = patch(Ex, Ey, 1, 'FaceColor', 'none');
end

