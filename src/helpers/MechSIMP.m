function [stiffness, stiffnessDer, thermalExp, thermalExpDer] = MechSIMP(materials, p_E, p_alpha)
%MECHSIMP Summary of this function goes here
%   Detailed explanation goes here
phi_min = 1e-12;
if numel(materials) < 2
    error("At least 2 materials must be specified");
elseif numel(materials) == 2
    stiffness = @(phi) materials(1).youngsModulus + phi.^p_E * (materials(2).youngsModulus - materials(1).youngsModulus);
    % Do not allow phi to be 0 when calculating the derivative
    % This to avoid Inf-values when p < 1
    stiffnessDer = @(phi) p_E * max(1e-12, phi).^(p_E-1) * (materials(2).youngsModulus - materials(1).youngsModulus);

    Dalpha = [materials.thermalExp] .* [materials.youngsModulus];

    thermalExp = @(phi) Dalpha(1) + phi.^p_alpha * (Dalpha(2) - Dalpha(1));
    thermalExpDer = @(phi) p_alpha * max(1e-12, phi).^(p_alpha-1) * (Dalpha(2) - Dalpha(1));
    
elseif numel(materials) == 3
    if numel(p_E) == 1
        p_E = p_E*ones(size(materials(2:end)));
    end
    if numel(p_alpha) == 1
        p_alpha = p_alpha*ones(size(materials(2:end)));
    end
    
    assert(numel(p_E) == 2 && numel(p_alpha) == 2, ...
        "The SIMP exponents must have the same dimensions as the desired design parameters");
    
    stiffness = @(phi) materials(1).youngsModulus + phi(1)^p_E(1)* ...
        (materials(2).youngsModulus + phi(2)^p_E(2)*(materials(3).youngsModulus - materials(2).youngsModulus) - ...
        materials(1).youngsModulus);
    % Do not allow phi to be 0 when calculating the derivative
    % This to avoid Inf-values when p < 1
    stiffnessDer = @(phi) [p_E(1)*max(phi_min, phi(1, :)).^(p_E(1)-1); 0] * ...
        (materials(2).youngsModulus - materials(1).youngsModulus) + ...
        [p_E(1)*max(phi_min, phi(1, :)).^(p_E(1)-1).*phi(2, :).^p_E(2); 
        phi(1, :).^p_E(2)*p_E(2).*max(phi_min, phi(2, :)).^(p_E(2)-1)] * ...
        (materials(3).youngsModulus - materials(2).youngsModulus);

    Dalpha = [materials.thermalExp] .* [materials.youngsModulus];

    thermalExp = @(phi) Dalpha(1) + phi(1, :).^p_alpha(1)* ...
        (Dalpha(2) + phi(2, :).^p_alpha(2)*(Dalpha(3) - Dalpha(2)) - ...
        Dalpha(1));
    thermalExpDer = @(phi) [p_alpha(1)*max(phi_min, phi(1, :)).^(p_alpha(1)-1); 0] * ...
        (Dalpha(2) - Dalpha(1)) + ...
        [p_alpha(1)*max(phi_min, phi(1, :)).^(p_alpha(1)-1)*phi(2, :).^p_alpha(2); 
        phi(1, :).^p_alpha(1)*p_alpha(2)*max(phi_min, phi(2, :))^(p_alpha(2)-1)] * ...
        (Dalpha(3) - Dalpha(2));
end
end

