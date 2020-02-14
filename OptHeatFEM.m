classdef OptHeatFEM < HeatFEM
    %HEATFEM Summary of this class goes here
    %   Detailed explanation goes here
    
    properties
        % Optimization properties
        mainEnod
        mainDensities
        mainVolumes
        
        conductivity
        heatCapacity
        density
    end
    
    methods
        function obj = OptHeatFEM(gmshData, tFinal, timeSteps)
            %HEATFEM Construct an instance of this class
            %   Detailed explanation goes here
            obj = obj@HeatFEM(gmshData, tFinal, timeSteps);
        end
        
        function addInterpFuncs(obj, conductivity, heatCapacity, density)
            % ADDINTERPFUNCS Add a interpolation function for the material
            % properties, with an element density as the only input
            obj.conductivity = conductivity;
            obj.heatCapacity = heatCapacity;
            obj.density = density;
        end
        
        function matrix = getElementBaseMatrix(obj, elementNbr, elementType, property)
            enod = obj.mainEnod(obj.mainEnod(:, 1) == elementNbr, 2:end);
            ex = obj.nodeCoordinates(enod, 1);
            ey = obj.nodeCoordinates(enod, 2);
            switch property
                case 'D'
                    matrix = obj.elementStiffness(ex, ey, elementType, ...
                        eye(2), obj.thickness);
                case 'cp'
                    matrix = obj.elementMass(ex, ey, elementType, 1, ...
                        obj.thickness);
                case 'rho'
                    matrix = obj.elementMass(ex, ey, elementType, 1, ...
                        obj.thickness);
            end
        end
        
        function [Ex, Ey, elementTemp] = getElemTemp(obj, timeStep)
            xs = obj.nodeCoordinates(:, 1);
            ys = obj.nodeCoordinates(:, 2);
            Ex = xs(obj.mainEnod(:, 2:end));
            Ey = ys(obj.mainEnod(:, 2:end));
            
            temp = obj.temperatures(:, timeStep+1);
            elementTemp = temp(obj.mainEnod(:, 2:end));
        end
        
        function reassemble(obj, densities)
            obj.mainDensities = densities;
            obj.reapplyBodyConditions();
            
            obj.loads = obj.fl + obj.fv + obj.fc;
            if obj.transient
                dt = obj.tFinal / (obj.timeSteps-1);
                obj.A = obj.M + dt*obj.theta*obj.K;
                obj.B = obj.M + dt*(obj.theta - 1)*obj.K;
                obj.loads = dt*obj.loads;
            end
        end
        
        function reapplyBodyConditions(obj)
            obj.K = sparse(obj.nbrDofs, obj.nbrDofs);
            obj.Kc = sparse(obj.nbrDofs, obj.nbrDofs);
            obj.M = sparse(obj.nbrDofs, obj.nbrDofs);
            for bc = obj.bodyConditions
                elementBlocks = getPhysicalEntity(obj.gmsh, bc.physicalName);
                switch bc.type
                    case 'main'
                        func = @(ex, ey, elementType, elementNbr) obj.elementStiffness(...
                            ex, ey, elementType, obj.elementProp('D', elementNbr), obj.thickness);
                        obj.K = obj.K + obj.integrate(elementBlocks, func, 2);
                        
                        func = @(ex, ey, elementType, elementNbr) obj.elementMass(...
                            ex, ey, elementType, obj.elementProp('rho', elementNbr)*...
                            obj.elementProp('cp', elementNbr), obj.thickness);
                        obj.M = obj.M + obj.integrate(elementBlocks, func, 2);
                end
            end
        end
    end
    methods(Access = protected)
        function [prop] = elementProp(obj, property, elementNbr)
            d = obj.mainDensities(obj.mainEnod(:, 1) == elementNbr);
            switch property
                case 'D'
                    prop = obj.conductivity(d);
                case 'cp'
                    prop = obj.heatCapacity(d);
                case 'rho'
                    prop = obj.density(d);
            end
        end
        
        function applyBodyConditions(obj)
            hasMain = 0;
            % Reset matrices
            obj.K = sparse(obj.nbrDofs, obj.nbrDofs);
            obj.Kc = sparse(obj.nbrDofs, obj.nbrDofs);
            obj.M = sparse(obj.nbrDofs, obj.nbrDofs);
            for bc = obj.bodyConditions
                elementBlocks = getPhysicalEntity(obj.gmsh, bc.physicalName);
                switch bc.type
                    case 'main'
                        hasMain = 1;
                        obj.mainEnod = getElements(obj.gmsh, elementBlocks(1).elementType);
                        obj.mainDensities = ones(size(obj.mainEnod, 1), 1);
                        obj.mainVolumes = obj.computeVolumes(obj.mainEnod);
                        
                        func = @(ex, ey, elementType, elementNbr) obj.elementStiffness(...
                            ex, ey, elementType, obj.elementProp('D', elementNbr), obj.thickness);
                        obj.K = obj.K + obj.integrate(elementBlocks, func, 2);
                        
                        func = @(ex, ey, elementType, elementNbr) obj.elementMass(...
                            ex, ey, elementType, obj.elementProp('rho', elementNbr)*...
                            obj.elementProp('cp', elementNbr), obj.thickness);
                        obj.M = obj.M + obj.integrate(elementBlocks, func, 2);
                    case 'load'
                        magnitude = bc.load;
                        func = @(ex, ey, elementType) obj.elementLoad(...
                            ex, ey, elementType, magnitude, obj.thickness);
                        f = obj.integrate(elementBlocks, func, 1);
                        obj.fv(:, bc.timeSteps) = ...
                            obj.fv(:, bc.timeSteps) + ...
                            repmat(f, [1 length(bc.timeSteps)]);
                end
            end
            if ~hasMain
                error('Model must have "main" body condition');
            end
        end
        
        function [matrix] = integrate(obj, elementBlocks, func, dim)
            if dim == 1
                matrix = zeros(obj.nbrDofs, 1);
            elseif dim == 2
                matrix = zeros(obj.nbrDofs, obj.nbrDofs);
            else
                error('The input "dim" must be either 1 or 2');
            end
            for block = elementBlocks
                for element = block.elements
                    ex = obj.nodeCoordinates(element.nodeTags, 1);
                    ey = obj.nodeCoordinates(element.nodeTags, 2);
                    edof = obj.Dofs(element.nodeTags);
                    if dim == 1
                        elemMatrix = func(ex, ey, block.elementType);
                        matrix(edof) = matrix(edof) + elemMatrix;
                    elseif dim == 2
                        elemMatrix = func(ex, ey, block.elementType, element.elementTag);
                        matrix(edof, edof) = matrix(edof, edof) + elemMatrix;
                    end
                end
            end
            matrix = sparse(matrix);
        end
        
        function volumes = computeVolumes(obj, Enod)
            volumes = zeros(size(Enod, 1), 1);
            for e = 1:size(Enod, 1)
                ex = obj.nodeCoordinates(Enod(e, 2:end), 1);
                ey = obj.nodeCoordinates(Enod(e, 2:end), 2);
                volumes(e) = polyarea(ex, ey) * obj.thickness;
            end
        end
    end
end
