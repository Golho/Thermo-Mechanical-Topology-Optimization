function [elementStresses] = smQUA_4stress(ex, ey, elementDisp, D, integrationRule, elementTemp, alpha)
%SMQUA_4STRESS Summary of this function goes here
%   Detailed explanation goes here
% ex - The element x coordinates, as a [4 x nbrElements] matrix
% ey - The element y coordinates, as a [4 x nbrElements] matrix
% elementDisp - The element displacements, as a [4 x nbrElements] matrix
% D - The coupling between the planar strains (xx, yy, xy) and a desired
% number of stresses
ex = reshape(ex, 4, []);
ey = reshape(ey, 4, []);
nbrElements = size(elementDisp, 2);
alphaExp = zeros(3, 1);
alphaExp(1:2) = alpha;

if nargin < 5
    integrationRule = 2;
end
if nargin < 6
    elementTemp = [];
end
nbrGaussPnts = integrationRule^2;
gaussPoints = zeros(nbrGaussPnts, 2);
weights = zeros(nbrGaussPnts, 1);

if integrationRule == 1
    oneDimGaussPnts = 0;
    oneDimWeights = 2;
elseif integrationRule == 2
    oneDimGaussPnts = [-0.577350269189626, 0.577350269189626];
    oneDimWeights = [1, 1];
elseif integrationRule == 3
    oneDimGaussPnts = [-0.7745966692414834, 0, 0.7745966692414834];
    oneDimWeights = [0.5555555555555556, 0.8888888888888888, 0.5555555555555556];
else
    error('Used number of integration points not implemented');
end

i = 1;
for x = oneDimGaussPnts
    for y = oneDimGaussPnts
        gaussPoints(i, :) = [x, y];
        i = i + 1;
    end
end

% Product of the weights in each dimension (for each point)
i = 1;
for wx = oneDimWeights
    for wy = oneDimWeights
        weights(i) = wx*wy;
        i = i + 1;
    end
end

xi = gaussPoints(:, 1);  
eta = gaussPoints(:, 2);

N = ones(nbrGaussPnts, 4);
% Assign to all gaussian points at the same time
N(:, [1 2]) = N(:, [1 2]) .* (1-eta);
N(:, [3 4]) = N(:, [3 4]) .* (1+eta);
N(:, [1 4]) = N(:, [1 4]) .* (1-xi);
N(:, [2 3]) = N(:, [2 3]) .* (1+xi);
N(:) = N / 4;

% Have every chunk of 3 rows below to one gaussian point
dNr = zeros(2*nbrGaussPnts, 4);
dNr(1:2:end, [1 3]) = dNr(1:2:end, [1 3]) + eta/4;
dNr(1:2:end, [2 4]) = dNr(1:2:end, [2 4]) - eta/4;
dNr(1:2:end, [2 3]) = dNr(1:2:end, [2 3]) + 1/4;
dNr(1:2:end, [1 4]) = dNr(1:2:end, [1 4]) - 1/4;

dNr(2:2:end, [1 3]) = dNr(2:2:end, [1 3]) + xi/4;
dNr(2:2:end, [2 4]) = dNr(2:2:end, [2 4]) - xi/4;
dNr(2:2:end, [1 2]) = dNr(2:2:end, [1 2]) - 1/4;
dNr(2:2:end, [3 4]) = dNr(2:2:end, [3 4]) + 1/4;

index = zeros(1, 2);
Be = zeros(3, 8);
elementStresses = zeros(size(D, 1), nbrGaussPnts*nbrElements);
Ne = zeros(1, 4);
for iElement = 1:nbrElements
    jacobianT = dNr*[ex(:, iElement), ey(:, iElement)];
    for i = 1:nbrGaussPnts
        index(:) = (1:2) + 2*(i-1);
        dNx = jacobianT(index, :) \ dNr(index, :);
        Be(1, 1:2:end) = dNx(1, :);
        Be(2, 2:2:end) = dNx(2, :);
        Be(3, 1:2:end) = dNx(2, :);
        Be(3, 2:2:end) = dNx(1, :);
        if ~isempty(elementTemp)
            Ne(:) = N(i, :);
            elementStrains = Be*elementDisp(:, iElement) - alphaExp*Ne*elementTemp(:, iElement);
        else
            elementStrains = Be*elementDisp(:, iElement);
        end
        
        elementStresses(:, nbrGaussPnts*(iElement-1) + i) = D*elementStrains;
    end
end
end

