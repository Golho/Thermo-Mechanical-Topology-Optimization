function [stiffness, stiffnessDer, thermalExp, thermalExpDer] = MechSIMP(material_1, material_2, p_E, p_alpha)
%MECHSIMP Summary of this function goes here
%   Detailed explanation goes here
stiffness = @(phi) material_1.youngsModulus + phi^p_E * (material_2.youngsModulus - material_1.youngsModulus);
% Do not allow phi to be 0 when calculating the derivative
% This to avoid Inf-values when p < 1
stiffnessDer = @(phi) p_E * max(1e-12, phi)^(p_E-1) * (material_2.youngsModulus - material_1.youngsModulus);

Dalpha_1 = material_1.thermalExp(1) * material_1.youngsModulus;
Dalpha_2 = material_2.thermalExp(1) * material_2.youngsModulus;

thermalExp = @(phi) Dalpha_1 + phi^p_alpha * (Dalpha_2 - Dalpha_1);
thermalExpDer = @(phi) p_alpha * max(1e-12, phi)^(p_alpha-1) * (Dalpha_2 - Dalpha_1);
end

