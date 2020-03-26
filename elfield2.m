function [p] = elfield2(Ex, Ey, Ed)
%ELFIELD2 Summary of this function goes here
%   Detailed explanation goes here
p = patch(Ex', Ey', Ed', 'EdgeColor', 'none');
axis("equal")
end

