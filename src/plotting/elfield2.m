function [p] = elfield2(Ex, Ey, Ed, varargin)
%ELFIELD2 Plot a scalar field on a 2D mesh
%   p = elfield2(Ex, Ey, Ed) Plots the elements defined by coordinates Ex
%   and Ey together with a scalar field Ed. The output is the patch object. 
%   Additional arguments to the patch function are carried by the varargin 
%   argument. 
p = patch(Ex, Ey, Ed, 'EdgeColor', 'none', varargin{:});
axis("equal");
end

