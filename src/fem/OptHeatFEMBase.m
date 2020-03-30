classdef OptHeatFEMBase < HeatFEMBase
    %   Detailed explanation goes here
    
    properties
        designPar
        
        conductivity
        heatCapacity
        
        conductivityDer
        heatCapacityDer
    end
    
    methods(Abstract)
        reassemble(obj, designPar)
        computeWeights(obj, radius)
        gradChainTerm(obj, adjointLoads, dT0dx)
    end
    
    methods
        function addInterpFuncs(obj, conductivity, conductivityDer, ...
                heatCapacity, heatCapacityDer)
            % ADDINTERPFUNCS Add a interpolation function for the material
            % properties, with an element design parameter as the only 
            % input
            obj.conductivity = conductivity;
            obj.heatCapacity = heatCapacity;
            obj.conductivityDer = conductivityDer;
            obj.heatCapacityDer = heatCapacityDer;
        end
        
        function matrix = getElementBaseMatrix(obj, ElementNbr, property)
            enod = obj.Enod(:, ElementNbr);
            elementCoord = obj.nodeCoordinates(:, enod)';
            switch property
                case 'D'
                    matrix = obj.elementStiffness(elementCoord, ...
                        obj.ElementType, eye(3));
                    
                case 'cp'
                    matrix = obj.elementMass(elementCoord, obj.ElementType, 1);

                case 'rho'
                    matrix = obj.elementMass(elementCoord, obj.ElementType, 1);
            end
        end
    end
    
    methods(Access = protected)
        function adjoints = solveAdjoint(obj, loads)
            if obj.transient
                % Create a matrix as big as the temperatures, but discard
                % the first column later (corresponding to lambda_N+1 which 
                % is not of interest
                adjoints = zeros(size(obj.temperatures));
                
                % Reverse the order of the loads and blocked dofs, to fit 
                % the transient solver
                loads = fliplr(loads);
                bds = fliplr(obj.blockedDofs);
                
                adjoints = obj.solveTransient(obj.A, obj.B, adjoints, ...
                    loads, bds);
                
                % Discard the first column corresponding to lambda_N+1 and
                % reverse the order of the adjoint column vectors
                adjoints = fliplr(adjoints(:, 2:end));
            else
                error('The steady state adjoint solver is not yet implemented');
            end
        end
        
        function [prop] = elementProp(obj, property, designPar)
            switch property
                case 'kappa'
                    if ~isempty(obj.conductivity)
                        prop = obj.conductivity(designPar);
                    else
                        prop = obj.D(1, 1);
                    end
                case 'D'
                    if ~isempty(obj.conductivity)
                        prop = obj.conductivity(designPar)*eye(3);
                    else
                        prop = obj.D;
                    end
                case 'cp'
                    if ~isempty(obj.heatCapacity)
                        prop = obj.heatCapacity(designPar);
                    else
                        prop = obj.cp;
                    end
                case 'rho'
                    if ~isempty(obj.density)
                        prop = obj.density(designPar);
                    else
                        prop = obj.rho;
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

