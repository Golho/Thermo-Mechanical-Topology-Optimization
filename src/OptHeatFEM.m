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
        
        function matrix = getElementBaseMatrix(obj, elementNbr, elementType, property, ex, ey)
            if nargin < 6
                enod = obj.mainEnod(obj.mainEnod(:, 1) == elementNbr, 2:end);
                ex = obj.nodeCoordinates(enod, 1);
                ey = obj.nodeCoordinates(enod, 2);
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
            % QUESTION: assign the variable "appliedBodies"?
            obj.K = sparse(obj.nbrDofs, obj.nbrDofs);
            obj.Kc = sparse(obj.nbrDofs, obj.nbrDofs);
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
            nbrMainElems = size(obj.mainEnod, 1);
            I = zeros(nbrMainElems*10, 1);
            J = zeros(nbrMainElems*10, 1);
            V = zeros(nbrMainElems*10, 1);
            elementCoord = zeros(nbrMainElems, 2);
            for e = 1:nbrMainElems
                enod = obj.mainEnod(e, 2:end);
                elementCoord(e, :) = mean(obj.nodeCoordinates(enod, 1:2), 1);
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
                
                adjoints = obj.solveTransient(adjoints, loads, bds);
                
                % Discard the first column corresponding to lambda_N+1 and
                % reverse the order of the adjoint column vectors
                adjoints = fliplr(adjoints(:, 2:end));
            else
                error('The steady state adjoint solver is not yet implemented');
            end
        end
    end
    methods(Access = protected)
        function [prop] = elementProp(obj, property, elementNbr, d)
            if nargin < 4
                d = obj.mainDensities(obj.mainEnod(:, 1) == elementNbr);
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
                case 'rho'
                    if ~isempty(obj.density)
                        prop = obj.density(d);
                    else
                        prop = obj.rho;
                    end
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
