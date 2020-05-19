function [Ke] = smQUA_4therm(ex, ey, D, alpha, thickness, integrationRule)
%SMQUA_4 Summary of this function goes here
%   Detailed explanation goes here
ex = reshape(ex, [], 1);
ey = reshape(ey, [], 1);
alpha = reshape(alpha, [], 1);

if nargin < 6
    integrationRule = 2;
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

% Have every chunk of 3 rows below to one gaussian point
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

Ke = zeros(8, 4);
JT = dNr*[ex, ey];
index = zeros(1, 2);
Be = zeros(3, 8);
Ne = zeros(1, 4);

for i = 1:nbrGaussPnts
    index(:) = (1:2) + 2*(i-1);
    detJ = abs(det(JT(index, :)));
    dNx = JT(index, :) \ dNr(index, :);
    Be(1, 1:2:end) = dNx(1, :);
    Be(2, 2:2:end) = dNx(2, :);
    Be(3, 1:2:end) = dNx(2, :);
    Be(3, 2:2:end) = dNx(1, :);
    Ne = N(i, :);
    Ke(:) = Ke + Be' * D(:, 1:2) * alpha * Ne * detJ * weights(i) * thickness;
end
end

