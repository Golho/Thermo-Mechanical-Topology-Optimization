function p = pathPlot(dofs, displacements, varargin)
%PATHPLOT Summary of this function goes here
%   Detailed explanation goes here
if length(dofs) == 2
    xDisp = displacements(dofs(1), :);
    yDisp = displacements(dofs(2), :);
    
    p = plot(xDisp, yDisp, varargin{:});
    axis equal
else
    warning("The plot for the specified spatial dimensino is not yet implemented");
end
end

