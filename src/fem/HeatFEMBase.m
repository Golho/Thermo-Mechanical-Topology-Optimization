classdef (Abstract) HeatFEMBase < FEMBase
    %UNTITLED Summary of this class goes here
    %   Detailed explanation goes here
    properties(Abstract)
       mesh 
    end
    
    properties
        temperatures
        nodeCoordinates
        fieldDim = 1;
        
        configuration
        
        M
        K
        theta
        
        % Debugging properties
        A % dR/dT_n matrix
        B % dR/dT_n-1 matrix
        loads
    end
    
    properties(Access = protected)
        fl
        fv
        blockedDofs
    end
    
    methods(Abstract)
        saveNodeField(obj, filePrefix, fieldVector, label)
        saveElementField(obj, filePrefix, fieldVector, label)
    end
    
    methods(Abstract, Access = protected)
        applyBoundaryConditions(obj)
        applyBodyConditions(obj)
    end
    
    methods
        function conf = get.configuration(obj)
            conf = obj.getConfiguration();
        end

        function assemble(obj)
            if ~obj.appliedBoundaries
                obj.applyBoundaryConditions();
            end
            if ~obj.appliedBodies
                obj.applyBodyConditions();
            end

            obj.loads = obj.fl + obj.fv;
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
                
                obj.temperatures = obj.solveTransient(obj.A, obj.B, ...
                    obj.temperatures, ...
                    weightedLoads, obj.blockedDofs);
            else
                obj.temperatures = obj.partitionAndSolve(obj.K, ...
                    obj.loads, obj.temperatures, obj.blockedDofs);
            end
        end
        
        function ed = getElemTemp(obj, timestep)
            temps = obj.temperatures(:, timestep+1);
            ed = temps(obj.Enod);
        end
    end
    
    methods(Access = protected)
        function conf = getConfiguration(obj)
            conf = struct(...
                'material', obj.material, ...
                'boundaryConditions', {obj.boundaryConditions}, ...
                'bodyConditions', {obj.bodyConditions}, ...
                'mesh', obj.mesh, ...
                'tFinal', obj.tFinal, ...
                'timeSteps', obj.timeSteps, ...
                'theta', obj.theta ...
                );
        end
        
        function init(obj)
            obj.Dofs = (1:obj.nbrDofs);
            % Initialize matrices
            obj.K = sparse(obj.nbrDofs, obj.nbrDofs);
            obj.M = sparse(obj.nbrDofs, obj.nbrDofs);
            
            % Initialize load vectors and temperature vector
            if obj.transient
                obj.fl = sparse(obj.nbrDofs, obj.timeSteps);
                obj.fv = sparse(obj.nbrDofs, obj.timeSteps);
                obj.loads = sparse(obj.nbrDofs, obj.timeSteps);
                obj.temperatures = zeros(obj.nbrDofs, obj.timeSteps);
                % Create cell array as each entry of blocked dofs may vary in size
                obj.blockedDofs = cell(1, obj.timeSteps);
            else
                obj.loads = sparse(obj.nbrDofs, 1);
                obj.temperatures = zeros(obj.nbrDofs, 1);
            end
        end
    end
    
    methods(Static)
        function [Ke] = elementStiffness(coord, elementType, D)
            switch elementType
                case Elements.TRI_3
                    Ke = hTRI_3K(coord(:, 1), coord(:, 2), D(1:2, 1:2), 1);
                case Elements.QUA_4 % QUADS
                    Ke = flw2i4e(coord(:, 1), coord(:, 2), [1 2], D(1:2, 1:2));
                    %Ke = hQUA_4K(ex, ey, obj.D, obj.thickness);
                case Elements.HEX_8
                    Ke = hHEX_8K(coord(:, 1), coord(:, 2), coord(:, 3), D, 2);
                otherwise
                    error("The element stiffness matrix is not yet implemented for the current element type");
            end
        end
        
        function [Me] = elementMass(coord, elementType, rho)
            switch elementType
                case Elements.TRI_3
                    Me = hTRI_3M(coord(:, 1), coord(:, 2), rho, 1);
                case Elements.QUA_4 % QUADS
                    Me = hQUA_4M(coord(:, 1), coord(:, 2), rho, 1);
                case Elements.HEX_8
                    Me = hHEX_8M(coord(:, 1), coord(:, 2), coord(:, 3), rho, 2);
                otherwise
                    error("The element mass matrix is not yet implemented for the current element type");
            end
        end
        
        function [fe] = elementLoad(coord, elementType, magnitude)
            switch elementType
                case Elements.PNT
                    fe = magnitude;
                case Elements.LIN_2
                    fe = hLIN_2f(coord(:, 1), coord(:, 2), magnitude, 1);
                otherwise
                    error("The element load vector is not yet implemented for the current element type");
            end
        end
        
        function X = solveTransient(A, B, X, F, bds)
            startTime = tic;
            % For each time step, solve the linear transient system AX=Y
            [nd, ts] = size(X);
            bd = reshape(bds{1}, [], 1);
            if isempty(bd)
                freeDofs = 1:nd;
            else
                freeDofs = ~any(1:nd == bd, 1);
            end
            % Decompose the A-matrix to speed up the linear system solver,
            % which prevents us from using the partitionAndSolve function
            % (Works only for linear systems, and for constant Dirichlet
            % boundary conditions)
            dA = decomposition(A(freeDofs, freeDofs));
            for timeStep = 1:(ts - 1)
                Yn = B * X(:, timeStep) + F(:, timeStep);
                Xn = X(:, timeStep+1);
                bd = reshape(bds{timeStep}, [], 1);
                
                if isempty(bd)
                    freeDofs = 1:nd;
                else
                    freeDofs = ~any(1:nd == bd, 1);
                end
                partYn = Yn(freeDofs) - ...
                    A(freeDofs, ~freeDofs)* ...
                    Xn(~freeDofs);
                
                partXn = dA \ partYn;
                Xn(freeDofs) = partXn;
                
                X(:, timeStep + 1) = Xn;
            end
            fprintf("Solved transient system:\t%f secs\n", toc(startTime));
        end
        
        function X = partitionAndSolve(A, Y, X, blockedDofs)
            % PARTITIONANDSOLVE Partition the system AX = Y and solve it
            bd = reshape(blockedDofs{1}, [], 1);
            if isempty(bd)
                freeDofs = 1:nd;
            else
                freeDofs = ~any(1:nd == bd, 1);
            end
            partY = Y(freeDofs) - ...
                A(freeDofs, ~freeDofs)* ...
                X(~freeDofs);
            partA = A(freeDofs, freeDofs);
            partX = partA \ partY;
            X(freeDofs) = partX;
        end
    end
end

