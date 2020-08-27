classdef Material
    %Material    Isotropic material with thermal and mechanical properties
    
    properties
        density
        heatCapacity
        conductivity        % thermal conductivity coefficient
        conductivityMatrix  % thermal conductivity matrix
        youngsModulus
        poissonsRatio
        stiffness           % material stiffness matrix
        thermalExp          % Thermal expansion coefficient
        thermalExpVector    % Thermal expansion vector
    end
    
    methods
        function obj = Material(density, heatCapacity, conductivity, E, nu, thermalExp)
            %Material Construct an isotropic material
            %   Constructor takes 6 scalar material properties:
            %       * Density
            %       * Heat Capacity
            %       * Conductivity
            %       * Young's Modulus
            %       * Poisson's ratio
            %       * Coefficient of Thermal Expansion (CTE)
            %   Properties left out automatically equals 1, except
            %   Poisson's ratio which will equal 0.25
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
            obj.youngsModulus = E;
            obj.poissonsRatio = nu;
            obj.thermalExp = thermalExp;
            
            % Create the material matrices from the isotropic material
            % properties
            obj.conductivityMatrix = conductivity*eye(3);
            obj.stiffness = isotropic(E, nu);
            obj.thermalExpVector = thermalExp*ones(3, 1);
        end
    end
end

