classdef (Abstract) OptMechFEMBase < MechFEMBase
    %OPTMECHFEMBASE Summary of this class goes here
    %   Detailed explanation goes here
    
    properties
        designPar
        
        stiffness
        thermalExp
        
        stiffnessDer
        thermalExpDer
        
        optConfiguration
    end
    
    methods(Abstract)
        reassemble(obj, designPar)
        computeWeights(obj, radius)
        gradChainTerm(obj, adjointLoads)
    end
    
    methods
        function conf = get.optConfiguration(obj)
            conf = struct( ...
                "stiffnessFunc", obj.stiffness, ...
                "thermalExpFunc", obj.thermalExp ...
            );
        end
        
        function addInterpFuncs(obj, stiffness, stiffnessDer, ...
                thermalExp, thermalExpDer)
            % ADDINTERPFUNCS Add a interpolation function for the material
            % properties, with an element design parameter as the only 
            % input
            obj.stiffness = stiffness;
            obj.thermalExp = thermalExp;
            obj.stiffnessDer = stiffnessDer;
            obj.thermalExpDer = thermalExpDer;
        end
        
        function matrix = getElementBaseMatrix(obj, ElementNbr, property)
            enod = obj.Enod(:, ElementNbr);
            elementCoord = obj.nodeCoordinates(:, enod)';
            switch property
                case 'D'
                    matrix = obj.elementStiffness(elementCoord, obj.ElementType, ...
                        isotropic(1, obj.material.poissonsRatio), obj.planarType);
                    
                case 'D-alpha'
                    matrix = obj.elementTempStiffness(elementCoord, obj.ElementType, ...
                        isotropic(1, obj.material.poissonsRatio), ones(3, 1), obj.planarType);
            end
        end
        
        function adjoints = solveAdjoint(obj, loads)
            % Create a matrix as big as the displacements
            adjoints = zeros(size(obj.displacements));
            
            adjoints = obj.partitionAndSolve(obj.K_tot, loads, adjoints, ...
                obj.blockedDofs);
        end
    end
    
    methods(Access = protected)
        function [prop] = elementProp(obj, property, designPar)
            switch property
                case 'E'
                    if ~isempty(obj.stiffness)
                        prop = obj.stiffness(designPar);
                    else
                        prop = obj.material.youngsModulus;
                    end
                case 'stiffness'
                    if ~isempty(obj.stiffness)
                        prop = obj.stiffness(designPar)*isotropic(1, obj.material.poissonsRatio);
                    else
                        prop = obj.material.stiffness;
                    end
                case 'alpha'
                    if ~isempty(obj.thermalExp)
                        prop = obj.thermalExp(designPar);
                    else
                        prop = obj.material.thermalExp;
                    end
                case 'alphavector'
                    if ~isempty(obj.thermalExp)
                        prop = obj.thermalExp(designPar)*ones(3, 1);
                    else
                        prop = obj.material.thermalExpVector;
                    end
            end
        end
    end
    
    methods(Static)
        function vol = elementVolume(ex, ey, ez, elementType)
            switch elementType
                case Elements.LIN_2
                    dx = ex(2, :) - ex(1, :);
                    dy = ey(2, :) - ey(1, :);
                    vol = sqrt(dx.^2 + dy.^2);
                case Elements.TRI_3
                    vol = polyarea(ex, ey);
                case Elements.QUA_4
                    vol = polyarea(ex, ey);
                case Elements.TET_4
                    vol = TET_4vol(ex, ey, ez);
                case Elements.HEX_8
                    vol = HEX_8vol(ex, ey, ez);
            end
        end
    end
end

