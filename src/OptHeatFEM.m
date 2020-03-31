classdef OptHeatFEM < HeatFEM
    %HEATFEM Summary of this class goes here
    %   Detailed explanation goes here
    
    properties
        % Optimization properties
        mainEnod
        mainDensities
        mainVolumes
        mainElementType
        
        conductivity
        heatCapacity
        
        conductivityDer
        heatCapacityDer
    end
    
    methods
        function obj = OptHeatFEM(gmshData, tFinal, timeSteps)
            %HEATFEM Construct an instance of this class
            %   Detailed explanation goes here
            obj = obj@HeatFEM(gmshData, tFinal, timeSteps);
        end
        
        function addInterpFuncs(obj, conductivity, conductivityDer, ...
                heatCapacity, heatCapacityDer)
            % ADDINTERPFUNCS Add a interpolation function for the material
            % properties, with an element density as the only input
            obj.conductivity = conductivity;
            obj.heatCapacity = heatCapacity;
            obj.conductivityDer = conductivityDer;
            obj.heatCapacityDer = heatCapacityDer;
        end
        
        function matrix = getMainElementBase(obj, mainElementNbr, property)
            enod = obj.mainEnod(2:end, mainElementNbr);
            ex = obj.nodeCoordinates(1, enod);
            ey = obj.nodeCoordinates(2, enod);
            matrix = obj.getElementBaseMatrix(0, obj.mainElementType, ...
                property, ex, ey);
        end
        
        function [Ex, Ey, elementTemp] = getMainElemTemp(obj, timeStep)
            xs = obj.nodeCoordinates(1, :);
            ys = obj.nodeCoordinates(2, :);
            Ex = xs(obj.mainEnod(2:end, :));
            Ey = ys(obj.mainEnod(2:end, :));
            
            temp = obj.temperatures(:, timeStep+1);
            elementTemp = temp(obj.mainEnod(2:end, :));
        end
        
        function reassemble(obj, densities)
            obj.mainDensities = densities;
            obj.reapplyBodyConditions();

            if obj.transient
                dt = obj.tFinal / (obj.timeSteps-1);
                obj.A = obj.M + dt*obj.theta*obj.K;
                obj.B = obj.M + dt*(obj.theta - 1)*obj.K;
            end
        end
        
        function reapplyBodyConditions(obj)
            obj.K = sparse(obj.nbrDofs, obj.nbrDofs);
            obj.M = sparse(obj.nbrDofs, obj.nbrDofs);
            for bc = obj.bodyConditions
                elementBlocks = getPhysicalEntity(obj.gmsh, bc.physicalName);
                switch bc.type
                    case 'main'
                        startTime = tic;
                        func = @(ex, ey, elementType, elementNbr) obj.elementProp('kappa', elementNbr)*...
                            obj.getElementBaseMatrix(elementNbr, elementType, 'D', ex, ey);
                        obj.K = obj.K + obj.integrate(elementBlocks, func, 2);
                        
                        func = @(ex, ey, elementType, elementNbr) obj.elementProp('cp', elementNbr)*...
                            obj.getElementBaseMatrix(elementNbr, elementType, 'cp', ex, ey);
                        obj.M = obj.M + obj.integrate(elementBlocks, func, 2);
                        fprintf("Reintegrated main BC:\t\t%f secs\n", toc(startTime));
                end
            end
        end
        
        function weights = computeMainWeights(obj, radius)
            startTime = tic;
            % Calculate the center point of each element
            nbrMainElems = size(obj.mainEnod, 2);
            I = zeros(nbrMainElems*10, 1);
            J = zeros(nbrMainElems*10, 1);
            V = zeros(nbrMainElems*10, 1);
            elementCoord = zeros(nbrMainElems, 2);
            for e = 1:nbrMainElems
                enod = obj.mainEnod(2:end, e);
                elementCoord(e, :) = mean(obj.nodeCoordinates(1:2, enod), 2);
            end
            % Calculate the weights
            c = 1;
            for e1 = 1:nbrMainElems
                c1 = elementCoord(e1, :);
                I(c) = e1;
                J(c) = e1;
                V(c) = radius;
                c = c + 1;
                for e2 = 1:(e1-1)
                    c2 = elementCoord(e2, :);
                    dx = c1(1) - c2(1);
                    dy = c1(2) - c2(2);
                    dist = sqrt(dx^2+dy^2);
                    if radius - dist > 0
                        I(c) = e1;
                        J(c) = e2;
                        V(c) = radius - dist;
                        c = c + 1;
                        I(c) = e2;
                        J(c) = e1;
                        V(c) = radius - dist;
                        c = c + 1;
                    end
                end
            end
            weights = sparse(I(1:c-1), J(1:c-1), V(1:c-1));
            
            % Normalize so the sum of the weights equal 1
            weights = weights ./ sum(weights, 2);
            fprintf("Computed weights:\t\t\t%f secs\n", toc(startTime));
        end
        
        function chainGrad = gradChainTerm(obj, adjointLoads, dT0dx)
            if nargin < 3
                dT0dx = sparse(size(obj.temperatures, 1), ...
                    length(obj.mainDensities));
            end
            
            deltaT = obj.tFinal / (obj.timeSteps - 1);
            
            adjoints = obj.solveAdjoint(adjointLoads);
            
            chainGrad = (-adjoints(:, 1)' * obj.B * dT0dx)';
            
            parfor e = 1:length(chainGrad)
                enod = obj.mainEnod(2:end, e);
                T_e = obj.temperatures(enod, :);
                adjoint_e = adjoints(enod, :);
                
                k0 = obj.getElementBaseMatrix(obj.mainEnod(1, e), ...
                    obj.mainElementType, 'D');
                c0 = obj.getElementBaseMatrix(obj.mainEnod(1, e), ...
                    obj.mainElementType, 'cp');
            
                dkappadphi = obj.conductivityDer(obj.mainDensities(e));
                dcpdphi = obj.heatCapacityDer(obj.mainDensities(e));

                dRdx = (deltaT*obj.theta*dkappadphi*k0 + dcpdphi*c0) * ...
                    T_e(:, 2:end) - ...
                    (deltaT*(obj.theta-1)*dkappadphi*k0 + dcpdphi*c0) * ...
                    T_e(:, 1:end-1);
                chainGrad(e) = chainGrad(e) - sum(dot(adjoint_e, dRdx));
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
                
                adjoints = obj.solveTransient(obj.A, obj.B, adjoints, loads, bds);
                
                % Discard the first column corresponding to lambda_N+1 and
                % reverse the order of the adjoint column vectors
                adjoints = fliplr(adjoints(:, 2:end));
            else
                error('The steady state adjoint solver is not yet implemented');
            end
        end
        
        function matrix = getElementBaseMatrix(obj, elementNbr, elementType, property, ex, ey)
            if nargin < 6
                enod = obj.mainEnod(2:end, obj.mainEnod(1, :) == elementNbr);
                ex = obj.nodeCoordinates(1, enod)';
                ey = obj.nodeCoordinates(2, enod)';
            end
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
        
        function [prop] = elementProp(obj, property, elementNbr, d)
            if nargin < 4
                d = obj.mainDensities(obj.mainEnod(1, :) == elementNbr);
            end
            switch property
                case 'kappa'
                    if ~isempty(obj.conductivity)
                        prop = obj.conductivity(d);
                    else
                        prop = obj.D(1, 1);
                    end
                case 'D'
                    if ~isempty(obj.conductivity)
                        prop = obj.conductivity(d)*eye(2);
                    else
                        prop = obj.D;
                    end
                case 'cp'
                    if ~isempty(obj.heatCapacity)
                        prop = obj.heatCapacity(d);
                    else
                        prop = obj.cp;
                    end
            end
        end
        
        function applyBodyConditions(obj)
            hasMain = false;
            % Reset matrices
            obj.K = sparse(obj.nbrDofs, obj.nbrDofs);
            obj.M = sparse(obj.nbrDofs, obj.nbrDofs);
            for bc = obj.bodyConditions
                elementBlocks = getPhysicalEntity(obj.gmsh, bc.physicalName);
                switch bc.type
                    case 'main'
                        hasMain = true;
                        obj.mainElementType = elementBlocks(1).elementType;
                        obj.mainEnod = getElements(obj.gmsh, elementBlocks(1).elementType);
                        obj.mainDensities = ones(size(obj.mainEnod, 2), 1);
                        obj.mainVolumes = obj.computeVolumes(obj.mainEnod);
                        
                        func = @(ex, ey, elementType, elementNbr) obj.elementStiffness(...
                            ex, ey, elementType, obj.elementProp('D', elementNbr), obj.thickness);
                        obj.K = obj.K + obj.integrate(elementBlocks, func, 2);
                        
                        func = @(ex, ey, elementType, elementNbr) obj.elementMass(...
                            ex, ey, elementType, obj.elementProp('cp', elementNbr), ...
                            obj.thickness);
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
        
        function volumes = computeVolumes(obj, Enod)
            volumes = zeros(size(Enod, 2), 1);
            for e = 1:size(Enod, 2)
                ex = obj.nodeCoordinates(1, Enod(2:end, e));
                ey = obj.nodeCoordinates(2, Enod(2:end, e));
                volumes(e) = polyarea(ex, ey) * obj.thickness;
            end
        end
    end
end
