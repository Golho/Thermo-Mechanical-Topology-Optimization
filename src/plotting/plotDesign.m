function [patchPlot] = plotDesign(Ex, Ey, designPar, patchPlot)
%PLOTDESIGN Summary of this function goes here
%   Detailed explanation goes here
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

