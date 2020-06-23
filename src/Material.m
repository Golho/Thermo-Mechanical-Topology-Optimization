classdef Material
    %UNTITLED2 Summary of this class goes here
    %   Detailed explanation goes here
    
    properties
        density
        heatCapacity
        conductivity % thermal conductivity coefficient
        conductivityMatrix % thermal conductivity matrix
        youngsModulus % Young's modulus
        poissonsRatio % Poisson's ratio
        stiffness % material stiffness matrix
        thermalExp % Thermal expansion coefficient
        thermalExpVector % Thermal expansion vector
    end
    
    methods
        function obj = Material(density, heatCapacity, conductivity, E, nu, thermalExp)
            %UNTITLED2 Construct an instance of this class
            %   Detailed explanation goes here
            if nargin < 6
                thermalExp = 1;
            end
            if nargin < 5
                E = 1;
            end
            if nargin < 4
                nu = 0.25;
            end
            if nargin < 3
                conductivity = 1;
            end
            if nargin < 2
                heatCapacity = 1;
            end
            if nargin < 1
                density = 1;
            end

            obj.density = density;
            obj.heatCapacity = heatCapacity;
            obj.conductivity = conductivity;
            obj.conductivityMatrix = conductivity*eye(3);
            obj.youngsModulus = E;
            obj.poissonsRatio = nu;
            obj.stiffness = isotropic(E, nu);
            obj.thermalExp = thermalExp;
            obj.thermalExpVector = thermalExp*ones(3, 1);
        end
    end
end

