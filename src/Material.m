classdef Material
    %UNTITLED2 Summary of this class goes here
    %   Detailed explanation goes here
    
    properties
        density
        heatCapacity
        Kappa % thermal conductivity matrix
        D % material stiffness matrix
    end
    
    methods
        function obj = Material(density, heatCapacity, Kappa, D)
            %UNTITLED2 Construct an instance of this class
            %   Detailed explanation goes here
            if nargin < 4
                D = eye(3);
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
            assert(size(D, 1) == 3 && size(D, 2), "The material stiffness matrix must be a 3x3 matrix");
            assert(size(Kappa, 1) == 3 && size(Kappa, 2), "The thermal conductivity matrix must be a 3x3 matrix");

            obj.density = density;
            obj.heatCapacity = heatCapacity;
            obj.Kappa = Kappa;
            obj.D = D;
        end
    end
end

