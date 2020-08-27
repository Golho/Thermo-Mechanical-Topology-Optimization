function [stiffnessMatrix] = isotropic(E, nu)
%ISOTROPIC Creates the material stiffness matrix for a linear isotropic
%material
%   stiffnessMatrix = isotropic(E, nu) Creates the stiffness matrix from
%   the Young's modulus E and Poisson's ratio nu
stiffnessMatrix = diag(ones(6, 1));
stiffnessMatrix(1:3, 1:3) = stiffnessMatrix(1:3, 1:3) + nu;
stiffnessMatrix = E / ((1 + nu)*(1-2*nu)) * (stiffnessMatrix - diag(2*nu*ones(6, 1)));
stiffnessMatrix(:, 4:6) = stiffnessMatrix(:, 4:6) / 2;
end

