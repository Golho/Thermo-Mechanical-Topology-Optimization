classdef HeatFEM < matlab.mixin.Copyable
    %HEATFEM Summary of this class goes here
    %   Detailed explanation goes here
    
    properties
        gmsh
        spatialDimensions
        nodeCoordinates
        Dofs
        % All the information about the mesh is given above
        % The properties below are just helpers
        nbrNodes
        nbrDofs
        
        % Boundary conditions
        nbrBcs = 0;
        boundaryConditions = struct(...
            'physicalName', [], ...
            'type', [], ...
            'value', [], ...
            'timeSteps', [] ...
        );
        
        % Body conditions
        nbrBodyCond = 0;
        bodyConditions = struct(...
            'type', [], ...
            'physicalName', [] ...
        );
    
        appliedBoundaries = 0;
        appliedBodies = 0;
        
        % Default material parameters
        thickness = 1;
        rho = 1;
        cp = 1;
        D = [1 0; 0 1]; % conductivity matrix
        
        A % dR/dT_n matrix
        B % dR/dT_n-1 matrix
        
        K % stiffness matrix
        M % mass matrix
        Kc % Robin stiffness matrix
        fc % Robin load vector
        fl % Neumann load vector
        fv % Volume load vector
        
        % A cell array of size [1 x N] where N is the number of time steps
        % Entry nbr n represents the blocked degrees of freedom at time t_n
        blockedDofs
        
        transient = 0;
        temperatures
        loads
        
        tFinal = -1;
        timeSteps = 0;
        theta = 1; % Time integration variable
    end
    
    methods
        function obj = HeatFEM(gmshData, tFinal, timeSteps)
            %HEATFEM Construct an instance of this class
            %   Detailed explanation goes here
            obj.gmsh = gmshData;
            obj.nodeCoordinates = getGlobalCoordinates(gmshData);
            
            obj.nbrNodes = size(obj.nodeCoordinates, 1);
            
            % TODO: Generalize for vector fields
            obj.nbrDofs = obj.nbrNodes;
            obj.Dofs = (1:obj.nbrNodes)';
            
            % Initialize matrices
            obj.K = sparse(obj.nbrDofs, obj.nbrDofs);
            obj.Kc = sparse(obj.nbrDofs, obj.nbrDofs);
            obj.M = sparse(obj.nbrDofs, obj.nbrDofs);
            
            if exist('tFinal', 'var')
                obj.transient = 1;
                obj.tFinal = tFinal;
                if ~exist('timeSteps', 'var')
                    error(['The number of time steps must be specified if ' ...
                        'a transient model is chosem']);
                end
                obj.timeSteps = timeSteps;
            end
            
            % Initialize load vectors and temperature vector
            if obj.transient
                obj.fl = zeros(obj.nbrDofs, obj.timeSteps);
                obj.fv = zeros(obj.nbrDofs, obj.timeSteps);
                obj.fc = zeros(obj.nbrDofs, obj.timeSteps);
                obj.loads = zeros(obj.nbrDofs, obj.timeSteps);
                obj.temperatures = zeros(obj.nbrDofs, obj.timeSteps);
                % Create cell array as each entry of blocked dofs may vary in size
                obj.blockedDofs = cell(1, timeSteps);
            else
                obj.loads = zeros(obj.nbrDofs, 1);
                obj.temperatures = zeros(obj.nbrDofs, 1);
            end
        end
        
        function setMaterial(obj, material, thickness)
            obj.D = material.D;
            obj.rho = material.density;
            obj.cp = material.heatCapacity;

            obj.thickness = thickness;
        end
        
        function addBoundaryCondition(obj, boundaryCondition)
            obj.boundaryConditions(obj.nbrBcs + 1) = boundaryCondition;
            obj.nbrBcs = obj.nbrBcs + 1;
        end
        
        function addBodyCondition(obj, bodyCondition)
            obj.bodyConditions(obj.nbrBodyCond + 1) = bodyCondition;
            obj.nbrBodyCond = obj.nbrBodyCond + 1;
        end
        
        function assemble(obj)
            obj.applyBoundaryConditions();
            obj.applyBodyConditions();
            
            obj.loads = obj.fl + obj.fv + obj.fc;
            if obj.transient
                dt = obj.tFinal / (obj.timeSteps-1);
                obj.A = obj.M + dt*obj.theta*obj.K;
                obj.B = obj.M + dt*(obj.theta - 1)*obj.K;
                obj.loads = dt*obj.loads;
            end
        end
        
        function solve(obj)
            if obj.transient
                % f_bar for n = 1, 2, ..., N
                weightedLoads = obj.theta*obj.loads(:, 2:end) + ...
                    (1-obj.theta)*obj.loads(:, 1:end-1);

                obj.temperatures = obj.solveTransient(obj.temperatures, ...
                    weightedLoads, obj.blockedDofs);
            else
                obj.temperatures = obj.partitionAndSolve(obj.K, ...
                    obj.loads, obj.temperatures, obj.blockedDofs);
            end
        end
        
        function [Ex, Ey, elementTemp] = getElemTemp(obj, elementType, timeStep)
            Enod = getElements(obj.gmsh, elementType);
            xs = obj.nodeCoordinates(:, 1);
            ys = obj.nodeCoordinates(:, 2);
            Ex = xs(Enod(:, 2:end));
            Ey = ys(Enod(:, 2:end));
            
            temp = obj.temperatures(:, timeStep+1);
            elementTemp = temp(Enod(:, 2:end));
        end
        
        function conf = getConfiguration(obj)
            conf = struct(...
                'thickness', obj.thickness, ...
                'rho', obj.rho, ...
                'cp', obj.cp, ...
                'D', obj.D, ...
                'boundaryConditions', obj.boundaryConditions, ...
                'bodyConditions', obj.bodyConditions, ...
                'mesh', obj.gmsh, ...
                'tFinal', obj.tFinal', ...
                'timeSteps', obj.timeSteps, ...
                'theta', obj.theta ...
            );
        end
    end
    
    methods(Access = protected)
        function applyBoundaryConditions(obj)
            if ~obj.appliedBoundaries
                for bc = obj.boundaryConditions
                    elementBlocks = getPhysicalEntity(obj.gmsh, bc.physicalName);
                    switch bc.type
                        case 'Neumann'
                            magnitude = bc.value;
                            func = @(ex, ey, elementType) obj.elementLoad(...
                                ex, ey, elementType, magnitude, obj.thickness);
                            f = obj.integrate(elementBlocks, func, 1);
                            obj.fl(:, bc.timeSteps) = ...
                                obj.fl(:, bc.timeSteps) + ...
                                repmat(f, [1 length(bc.timeSteps)]);
                        case 'Robin'
                            magnitude = bc.alpha * bc.Tinf;
                            func = @(ex, ey, elementType) obj.elementLoad(...
                                ex, ey, elementType, magnitude, obj.thickness);
                            f = obj.integrate(elementBlocks, func, 1);
                            obj.fc(:, bc.timeSteps) = ...
                                obj.fc(:, bc.timeSteps) + ...
                                repmat(f, [1 length(bc.timeSteps)]);

                            magnitude = bc.alpha;
                            func = @(ex, ey, elementType) obj.elementMass(...
                                ex, ey, elementType, magnitude, obj.thickness);
                            % TODO: handle different Kc at different time steps
                            obj.Kc = obj.Kc + obj.integrate(elementBlocks, func, 2);

                        case 'Dirichlet'
                            % Extract all the nodes of the element blocks
                            nodes = [];
                            for elementBlock = elementBlocks
                                nodes = horzcat(nodes, [elementBlock.elements.nodeTags]);
                            end
                            bd = obj.Dofs(unique(nodes));
                            addBd = @(c) {[c, bd]};
                            obj.blockedDofs(bc.timeSteps) = cellfun(addBd, ...
                                obj.blockedDofs(bc.timeSteps));
                            obj.temperatures(bd, bc.timeSteps) = bc.value;
                    end
                end
                obj.appliedBoundaries = 1;
            end
        end
        
        function applyBodyConditions(obj)
            if ~obj.appliedBodies
                for bc = obj.bodyConditions
                    elementBlocks = getPhysicalEntity(obj.gmsh, bc.physicalName);
                    switch bc.type
                        case 'main'
                            func = @(ex, ey, elementType) obj.elementStiffness(...
                                ex, ey, elementType, obj.D, obj.thickness);
                            obj.K = obj.K + obj.integrate(elementBlocks, func, 2);

                            magnitude = obj.rho * obj.cp;
                            func = @(ex, ey, elementType) obj.elementMass(...
                                ex, ey, elementType, magnitude, obj.thickness);
                            obj.M = obj.M + obj.integrate(elementBlocks, func, 2);
                        case 'load'
                            magnitude = bc.value;
                            func = @(ex, ey, elementType) obj.elementLoad(...
                                ex, ey, elementType, magnitude, obj.thickness);
                            f = obj.integrate(elementBlocks, func, 1);
                            obj.fv(:, bc.timeSteps) = ...
                                obj.fv(:, bc.timeSteps) + ...
                                repmat(f, [1 length(bc.timeSteps)]);
                        case 'initial'
                            % Extract all the nodes of the element blocks
                            nodes = [];
                            for elementBlock = elementBlocks
                                nodes = horzcat(nodes, [elementBlock.elements.nodeTags]);
                            end
                            bd = obj.Dofs(unique(nodes));
                            obj.temperatures(bd, 1) = bc.value;
                    end
                end
                obj.appliedBodies = 1;
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
                    elemMatrix = func(ex, ey, block.elementType);
                    if dim == 1
                        matrix(edof) = matrix(edof) + elemMatrix;
                    elseif dim == 2
                        matrix(edof, edof) = matrix(edof, edof) + elemMatrix;
                    end
                end
            end
            matrix = sparse(matrix);
        end
        
        function X = solveTransient(obj, X, F, bds)
            % For each time step, solve the linear transient system AX=Y
            bd = bds{1};
            if isempty(bd)
                freeDofs = 1:obj.nbrDofs;
            else
                freeDofs = ~any(1:obj.nbrDofs == bd);
            end
            A = sparse(obj.A);
            % Decompose the A-matrix to speed up the linear system solver,
            % which prevents us from using the partitionAndSolve function
            % (Works only for linear systems, and for constant Dirichlet
            % boundary conditions)
            dA = decomposition(A(freeDofs, freeDofs));
            for timeStep = 1:(obj.timeSteps - 1)
                Yn = obj.B * X(:, timeStep) + F(:, timeStep);
                Xn = X(:, timeStep+1);
                bd = bds{timeStep};
                
                if isempty(bd)
                    freeDofs = 1:obj.nbrDofs;
                else
                    freeDofs = ~any(1:obj.nbrDofs == bd);
                end
                partYn = Yn(freeDofs) - ...
                        A(freeDofs, ~freeDofs)* ...
                        Xn(~freeDofs);
                
                partXn = dA \ partYn;
                Xn(freeDofs) = partXn;
                
                X(:, timeStep + 1) = Xn;
            end
        end
    end
    
    methods(Static)
        function [Ke] = elementStiffness(ex, ey, elementType, D, thickness)
            switch elementType
                case Elements.TRI_3
                    Ke = hTRI_3K(ex, ey, D, thickness);
                case Elements.QUA_4 % QUADS
                    Ke = flw2i4e(ex', ey', [thickness 2], D);
                    %Ke = hQUA_4K(ex, ey, obj.D, obj.thickness);
                otherwise
                    error("The element stiffness matrix is not yet implemented for the current element type");
            end
        end
        
        function [Me] = elementMass(ex, ey, elementType, rho, thickness)
            switch elementType
                case Elements.TRI_3
                    Me = hTRI_3M(ex, ey, rho, thickness);
                case Elements.QUA_4 % QUADS
                    Me = hQUA_4M(ex, ey, rho, thickness);
                otherwise
                    error("The element mass matrix is not yet implemented for the current element type");
            end
        end
        
        function [fe] = elementLoad(ex, ey, elementType, magnitude, thickness)
            switch elementType
                case Elements.PNT
                    fe = magnitude*thickness;
                case Elements.LIN_2
                    fe = hLIN_2f(ex, ey, magnitude, thickness);
                otherwise
                    error("The element load vector is not yet implemented for the current element type");
            end
        end
        
        function [X] = partitionAndSolve(A, Y, X, blockedDofs)
            % PARTITIONANDSOLVE Partition the system AX = Y and solve it
            freeDofs = ~any(1:length(X) == blockedDofs);
            partY = Y(freeDofs) - ...
                    A(freeDofs, ~freeDofs)* ...
                    X(~freeDofs);
            partA = A(freeDofs, freeDofs);
            partX = partA \ partY;
            X(freeDofs) = partX;
        end
    end
end
