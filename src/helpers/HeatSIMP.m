function [conductivity, conductivityDer, heatCapacity, heatCapacityDer] = HeatSIMP(materials, p_kappa, p_cp)
%MECHSIMP Summary of this function goes here
%   Detailed explanation goes here
phi_min = 1e-12;
if numel(materials) < 2
    error("At least 2 materials must be specified");
elseif numel(materials) == 2
    conductivity = @(phi) materials(1).conductivity + phi.^p_kappa * (materials(2).conductivity - materials(1).conductivity);
    % Do not allow phi to be 0 when calculating the derivative
    % This to avoid Inf-values when p < 1
    conductivityDer = @(phi) p_kappa * max(1e-12, phi).^(p_kappa-1) * (materials(2).conductivity - materials(1).conductivity);

    volumetricHeat = [materials.density] .* [materials.heatCapacity];

    heatCapacity = @(phi) volumetricHeat(1) + phi.^p_cp * (volumetricHeat(2) - volumetricHeat(1));
    heatCapacityDer = @(phi) p_cp * max(1e-12, phi).^(p_cp-1) * (volumetricHeat(2) - volumetricHeat(1));
    
elseif numel(materials) == 3
    if numel(p_kappa) == 1
        p_kappa = p_kappa*ones(size(materials(2:end)));
    end
    if numel(p_cp) == 1
        p_cp = p_cp*ones(size(materials(2:end)));
    end
    assert(numel(p_kappa) == 2 && numel(p_cp) == 2, ...
        "The SIMP exponents must have the same dimensions as the desired design parameters");
    
    conductivity = @(phi) materials(1).conductivity + phi(1, :).^p_kappa(1).* ...
        (materials(2).conductivity + phi(2, :).^p_kappa(2)*(materials(3).conductivity - materials(2).conductivity) - ...
        materials(1).conductivity);
    % Do not allow phi to be 0 when calculating the derivative
    % This to avoid Inf-values when p < 1
    conductivityDer = @(phi) [p_kappa(1)*max(phi_min, phi(1, :)).^(p_kappa(1)-1); 0] * ...
        (materials(2).conductivity - materials(1).conductivity) + ...
        [p_kappa(1)*max(phi_min, phi(1, :)).^(p_kappa(1)-1).*phi(2, :).^p_kappa(2); 
        phi(1, :).^p_kappa(2)*p_kappa(2).*max(phi_min, phi(2, :)).^(p_kappa(2)-1)] * ...
        (materials(3).conductivity - materials(2).conductivity);

    volumetricHeat = [materials.density] .* [materials.heatCapacity];

    heatCapacity = @(phi) volumetricHeat(1) + phi(1, :).^p_cp(1)* ...
        (volumetricHeat(2) + phi(2, :).^p_cp(2)*(volumetricHeat(3) - volumetricHeat(2)) - ...
        volumetricHeat(1));
    heatCapacityDer = @(phi) [p_cp(1)*max(phi_min, phi(1, :)).^(p_cp(1)-1); 0] * ...
        (volumetricHeat(2) - volumetricHeat(1)) + ...
        [p_cp(1)*max(phi_min, phi(1, :)).^(p_cp(1)-1)*phi(2, :).^p_cp(2); 
        phi(1, :).^p_cp(1)*p_cp(2)*max(phi_min, phi(2, :))^(p_cp(2)-1)] * ...
        (volumetricHeat(3) - volumetricHeat(2));
end
end

