function [plotObj1, plotObj2] = plotIntermediate(femModel, designPar, figHandle1, figHandle2, plotObj1, plotObj2)
%PLOTINTERMEDIATE Summary of this function goes here
%   Detailed explanation goes here
[Ex, Ey, ed] = femModel.getElemTemp(femModel.timeSteps-1);
if nargin < 5
    figure(figHandle1)
    plotObj1 = elfield2(Ex, Ey, ed);
    figure(figHandle2)
    plotObj2 = elfield2(Ex, Ey, designPar);
else
    plotObj1.CData = ed';
    if nargin >= 6
        plotObj2.CData = designPar;
    end
    drawnow
end
end

