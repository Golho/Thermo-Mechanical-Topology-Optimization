function [p] = eldisplace2(Ex, Ey, ed, scaleFactor, varargin)
%ELDISPLACE2 Plot a deformed 2D mesh
%   plot = eldisplace2(Ex, Ey, ed, scaleFactor) Plot the elements with
%   element coordinates Ex, Ey and displacements ed. The displacements are
%   scaled with scaleFactor. The output is the patch object. Additional
%   arguments to the patch function are carried by the varargin argument.
Ex = Ex + scaleFactor*ed(1:2:end, :);
Ey = Ey + scaleFactor*ed(2:2:end, :);

p = patch(Ex, Ey, 1, 'FaceColor', 'none', varargin{:});
end

