function [conductivity, conductivityDer, heatCapacity, heatCapacityDer] = HeatSIMP(material_1, material_2, p_kappa, p_cp)
%MECHSIMP Summary of this function goes here
%   Detailed explanation goes here
conductivity = @(phi) material_1.Kappa(1) + phi^p_kappa * (material_2.Kappa(1)- material_1.Kappa(1));
conductivityDer = @(phi) p_kappa * phi^(p_kappa-1) * (material_2.Kappa(1) - material_1.Kappa(1));

heatCapacity = @(phi) material_1.heatCapacity + phi^p_cp * (material_2.heatCapacity - material_1.heatCapacity);
heatCapacityDer = @(phi) p_cp * phi^(p_cp-1) * (material_2.heatCapacity - material_1.heatCapacity);
end

