function [patchPlot] = plotDesign(Ex, Ey, designPar, patchPlot)
%PLOTDESIGN Plot a topology optimized design
%   patchPlot = plotDesign(Ex, Ey, designPar) Plot the design
%   parameters in elements defined by element coordinates Ex and Ey. The
%   design is defined in designPar as floats between 0 and 1. The function
%   accepts either 1 or 2 design parameters per element. 
%   
%   patchPlot = plotDesign(Ex, Ey, designPar, patchPlot) Update the
%   patch plot patchPlot
if nargin < 4
    if size(designPar, 1) == 2
        patchPlot = elfield2(Ex, Ey, designPar(2, :));
        patchPlot.FaceAlpha = "flat";
        patchPlot.FaceVertexAlphaData = designPar(1, :)';
        patchPlot.AlphaDataMapping = "none";
    elseif size(designPar, 1) == 1
        patchPlot = elfield2(Ex, Ey, designPar);
    end
    colorbar;
    caxis([0, 1]);
else
    % Only updates patchPlot
    if size(designPar, 1) == 2
        patchPlot.CData = designPar(2, :);
        patchPlot.FaceVertexAlphaData = designPar(1, :)';
    elseif size(designPar, 1) == 1
        patchPlot.CData = designPar;
    end
end
drawnow;
end

