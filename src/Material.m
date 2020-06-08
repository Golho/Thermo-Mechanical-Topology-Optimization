classdef Material
    %UNTITLED2 Summary of this class goes here
    %   Detailed explanation goes here
    
    properties
        density
        heatCapacity
        Kappa % thermal conductivity matrix
        youngsModulus % Young's modulus
        poissonsRatio % Poisson's ratio
        stiffness % material stiffness matrix
        thermalExp % Thermal expansion coefficient
        thermalExpVector % Thermal expansion vector
    end
    
    methods
        function obj = Material(density, heatCapacity, Kappa, E, nu, thermalExp)
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
                Kappa = eye(3);
            end
            if nargin < 2
                heatCapacity = 1;
            end
            if nargin < 1
                density = 1;
            end
            assert(size(Kappa, 1) == 3 && size(Kappa, 2), "The thermal conductivity matrix must be a 3x3 matrix");

            obj.density = density;
            obj.heatCapacity = heatCapacity;
            obj.Kappa = Kappa;
            obj.youngsModulus = E;
            obj.poissonsRatio = nu;
            obj.stiffness = isotropic(E, nu);
            obj.thermalExp = thermalExp;
            obj.thermalExpVector = thermalExp*ones(3, 1);
        end
    end
end

