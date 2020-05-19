function [stiffnessMatrix] = isotropic(E, nu)
%UNTITLED2 Summary of this function goes here
%   Detailed explanation goes here
stiffnessMatrix = diag(ones(6, 1));
stiffnessMatrix(1:3, 1:3) = stiffnessMatrix(1:3, 1:3) + nu;
stiffnessMatrix = E / ((1 + nu)*(1-2*nu)) * (stiffnessMatrix - diag(2*nu*ones(6, 1)));
stiffnessMatrix(:, 4:6) = stiffnessMatrix(:, 4:6) / 2;
end

