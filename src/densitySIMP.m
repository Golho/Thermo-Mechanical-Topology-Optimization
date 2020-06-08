function [density, densityDer] = densitySIMP(materials, p)
%DENSITYSIMP Summary of this function goes here
%   Detailed explanation goes here
phi_min = 1e-12;
if numel(materials) < 2
    error("At least 2 materials must be specified");
elseif numel(materials) == 2
    density = @(phi) materials(1).density + phi.^p * (materials(2).density - materials(1).density);
    % Do not allow phi to be 0 when calculating the derivative
    % This to avoid Inf-values when p < 1
    densityDer = @(phi) p * max(1e-12, phi).^(p-1) * (materials(2).density - materials(1).density);
    
elseif numel(materials) == 3
    if numel(p) < 2
        p = p*ones(2, 1);
    end
    
    density = @(phi) materials(1).density + phi(1, :).^p(1).* ...
        (materials(2).density + phi(2, :).^p(2)*(materials(3).density - materials(2).density) - ...
        materials(1).density);
    % Do not allow phi to be 0 when calculating the derivative
    % This to avoid Inf-values when p < 1
    densityDer = @(phi) [p(1)*max(phi_min, phi(1, :)).^(p(1)-1); zeros(1, size(phi, 2))] * ...
        (materials(2).density - materials(1).density) + ...
        [p(1)*max(phi_min, phi(1, :)).^(p(1)-1).*phi(2, :).^p(2); 
        phi(1, :).^p(2)*p(2).*max(phi_min, phi(2, :)).^(p(2)-1)] * ...
        (materials(3).density - materials(2).density);
end
end

