function [Ke] = hHEX_8K(ex, ey, ez, D, integrationRule)
ex = reshape(ex, [], 1);
ey = reshape(ey, [], 1);
ez = reshape(ez, [], 1);

if nargin < 5
    integrationRule = 2;
end
nbrGaussPnts = integrationRule^3;
gaussPoints = zeros(nbrGaussPnts, 3);
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
        for z = oneDimGaussPnts
            gaussPoints(i, :) = [x, y, z];
            i = i + 1;
        end
    end
end

% Product of the weights in each dimension (for each point)
i = 1;
for wx = oneDimWeights
    for wy = oneDimWeights
        for wz = oneDimWeights
            weights(i) = wx*wy*wz;
            i = i + 1;
        end
    end
end

xi = gaussPoints(:,1);
eta = gaussPoints(:,2);
zeta = gaussPoints(:, 3);

% Have every chunk of 3 rows below to one gaussian point
dNr = zeros(3*nbrGaussPnts, 8);
s(:, 1) = (xi-1).*(zeta-1) / 8;
s(:, 2) = (xi-1).*(zeta+1) / 8;
s(:, 3) = (xi+1).*(zeta-1) / 8;
s(:, 4) = (xi+1).*(zeta+1) / 8;
s(:, 5) = (eta-1).*(zeta-1) / 8;
s(:, 6) = (eta-1).*(zeta+1) / 8;
s(:, 7) = (eta+1).*(zeta-1) / 8;
s(:, 8) = (eta+1).*(zeta+1) / 8;
s(:, 9) = (eta-1).*(xi-1) / 8;
s(:, 10) = (eta-1).*(xi+1) / 8;
s(:, 11) = (eta+1).*(xi-1) / 8;
s(:, 12) = (eta+1).*(xi+1) / 8;
% Assign to all gaussian points at the same time (thereof 1:3:end)
dNr(1:3:end, :) = [-s(:, 5), s(:, 5), -s(:, 7), s(:, 7), ...
                   s(:, 6), -s(:, 6), s(:, 8), -s(:, 8)];

dNr(2:3:end, :) = [-s(:, 1), s(:, 3), -s(:, 3), s(:, 1), ...
                   s(:, 2), -s(:, 4), s(:, 4), -s(:, 2)];

dNr(3:3:end, :) = [-s(:, 9), s(:,10), -s(:,12), s(:,11), ...
                   s(:, 9), -s(:,10), s(:,12), -s(:,11)];

Ke = zeros(8,8);
B = zeros(3, 8);
JT = dNr*[ex, ey, ez];
index = zeros(1, 3);

for i = 1:nbrGaussPnts
    index(:) = (1:3) + 3*(i-1);
    detJ = abs(det(JT(index, :)));
    B(:) = JT(index, :) \ dNr(index, :);
    Ke(:) = Ke + B' * D * B * detJ * weights(i);
end
end
