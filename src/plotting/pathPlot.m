function p = pathPlot(dofs, displacements, varargin)
%PATHPLOT Plot the path from displacements
%   p = patchPlot(dofs, displacments) Plot the path defined in
%   displacements, only defined by a two-element vector dofs containing
%   degrees of freedom in x and y to plot.
if length(dofs) == 2
    xDisp = displacements(dofs(1), :);
    yDisp = displacements(dofs(2), :);
    
    p = plot(xDisp, yDisp, varargin{:});
    axis equal
else
    warning("The plot for the specified spatial dimension is not yet implemented");
end
end

