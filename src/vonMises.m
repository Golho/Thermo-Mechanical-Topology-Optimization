function [vmStresses] = vonMises(stresses)
%VONMISES Summary of this function goes here
%   Detailed explanation goes here
numStresses = size(stresses, 1);
if numStresses == 3
    vmStresses = sqrt(stresses(1, :).^2 - stresses(1, :).*stresses(2, :) + ...
        stresses(2, :).^2 + 3*stresses(3, :).^2);
end
end

