classdef (Abstract) MechFEMBase < FEMBase
    %UNTITLED2 Summary of this class goes here
    %   Detailed explanation goes here
    
    properties(Abstract, Access = protected)
        mesh
    end
    
    properties
        temperatureChanges
        displacements
        nodeCoordinates
        fieldDim;
        
        configuration
        
        K
        loads
        K_thermal
        % Debugging properties
        
    end
    
    properties(Access = protected)
        K_c
        fl % boundary load
        fv % volume load
        ft % thermal load
        blockedDofs
        planarType % Plane stress or plane strain if 2D
    end
    
    properties(Dependent)
        K_tot
    end
    
    methods(Abstract)
        saveNodeField(obj, filePrefix, fieldVector, label)
        saveNodeVectorField(obj, filePrefix, fieldVector, label)
        saveElementField(obj, filePrefix, fieldVector, label)
        f = getDummy(obj, name)
    end
    
    methods(Abstract, Access = protected)
        applyBoundaryConditions(obj)
        applyBodyConditions(obj)
    end
    
    methods
        function K_tot = get.K_tot(obj)
            K_tot = obj.K + obj.K_c;
        end
        
        function conf = get.configuration(obj)
            conf = struct(...
                'material', obj.material, ...
                'boundaryConditions', {obj.boundaryConditions}, ...
                'bodyConditions', {obj.bodyConditions}, ...
                'mesh', obj.mesh, ...
                'timeSteps', obj.timeSteps ...
                );
        end
        
        function assemble(obj)
            if ~obj.appliedBoundaries
                obj.applyBoundaryConditions();
            end
            if ~obj.appliedBodies
                obj.applyBodyConditions();
            end
            
            obj.loads = obj.fl + obj.fv + obj.ft;
        end
        
        function setTemperatures(obj, temperatureChanges)
            obj.temperatureChanges = temperatureChanges;
            obj.ft = obj.K_thermal * obj.temperatureChanges;
            obj.loads = obj.fl + obj.fv + obj.ft;
        end
        
        function solve(obj)
            obj.displacements = obj.partitionAndSolve(obj.K_tot, ...
                obj.loads, obj.displacements, obj.blockedDofs);
        end
        
        function ed = getElemDisp(obj, timeStep)
            timeDisp = obj.displacements(:, timeStep+1);
            ed = timeDisp(obj.Edof);
        end
        
        function timeStress = getElemStress(obj, timeStep)
            timeDisp = obj.getElemDisp(timeStep);
            if isempty(obj.temperatureChanges)
                timeStress = obj.elementStress(timeDisp);
            else
                et = obj.temperatureChanges(:, timeStep+1);
                timeTemp = et(obj.Enod);
                timeStress = obj.elementStress(timeDisp, timeTemp);
            end
        end
    end
    
    methods(Access = protected)
        function init(obj)
            obj.Dofs = zeros(obj.fieldDim, obj.nbrNodes);
            obj.Dofs(:) = 1:obj.nbrDofs;
            % Initialize matrices
            obj.K = sparse(obj.nbrDofs, obj.nbrDofs);
            obj.K_c = sparse(obj.nbrDofs, obj.nbrDofs);
            
            obj.fl = sparse(obj.nbrDofs, obj.timeSteps);
            obj.fv = sparse(obj.nbrDofs, obj.timeSteps);
            obj.ft = sparse(obj.nbrDofs, obj.timeSteps);
            obj.loads = sparse(obj.nbrDofs, obj.timeSteps);
            obj.displacements = zeros(obj.nbrDofs, obj.timeSteps);
            % Create cell array as each entry of blocked dofs may vary in size
            obj.blockedDofs = cell(1, obj.timeSteps);
        end
        
        function [es] = elementStress(obj, elementDisp, elementTemp)
            switch obj.ElementType
                case Elements.QUA_4 % QUADS
                    switch obj.planarType
                        case "plane strain"
                            D = obj.material.stiffness([1, 2, 4], [1, 2, 4]);
                        case "plane stress"
                            C = inv(obj.material.stiffness);
                            D = inv(C([1 2 4], [1 2 4]));
                        otherwise
                            error("Planar type must be given for planar elements")
                    end
                    if exist('elementTemp', 'var')
                        es = smQUA_4stress(obj.Ex, obj.Ey, elementDisp, D, ...
                            1, elementTemp, obj.material.thermalExp(1:2));
                    else
                        es = smQUA_4stress(obj.Ex, obj.Ey, elementDisp, D, 1);
                    end
                    %Ke = hQUA_4K(ex, ey, obj.D, obj.thickness);
                otherwise
                    error("The element stiffness matrix is not yet implemented for the current element type");
            end
        end
    end
    
    methods(Static)
        function [Ke] = elementStiffness(coord, elementType, stiffness, planarType)
            if nargin < 4
                planarType = "plane strain";
            end
            switch elementType
                case Elements.QUA_4 % QUADS
                    switch planarType
                        case "plane strain"
                            DD = stiffness([1, 2, 4], [1, 2, 4]);
                        case "plane stress"
                            C = inv(stiffness);
                            DD = inv(C([1 2 4], [1 2 4]));
                        otherwise
                            error("Planar type must be given for planar elements")
                    end
                    Ke = smQUA_4K(coord(:, 1), coord(:, 2), DD, 1);
                    %Ke = hQUA_4K(ex, ey, obj.D, obj.thickness);
                otherwise
                    error("The element stiffness matrix is not yet implemented for the current element type");
            end
        end
        
        function [Ke] = elementTempStiffness(coord, elementType, stiffness, thermalExp, planarType)
            if nargin < 5
                planarType = "plane strain";
            end
            switch elementType
                case Elements.QUA_4 % QUADS
                    switch planarType
                        case "plane strain"
                            D = stiffness([1, 2, 4], [1, 2, 4]);
                        case "plane stress"
                            C = inv(stiffness);
                            D = inv(C([1 2 4], [1 2 4]));
                        otherwise
                            error("Planar type must be given for planar elements")
                    end
                    Ke = smQUA_4therm(coord(:, 1), coord(:, 2), D, thermalExp(1:2), 1);
                    %Ke = hQUA_4K(ex, ey, obj.D, obj.thickness);
                otherwise
                    error("The element stiffness matrix is not yet implemented for the current element type");
            end
        end
        
        function [fe] = elementLoad(coord, elementType, magnitudes)
            switch elementType
                case Elements.PNT
                    fe = magnitudes;
                case Elements.LIN_2
                    fe = smLIN_2f(coord(:, 1), coord(:, 2), magnitudes, 1);
                otherwise
                    error("The element load vector is not yet implemented for the current element type");
            end
        end
        
        function X = partitionAndSolve(A, Y, X, blockedDofs)
            % PARTITIONANDSOLVE Partition the system AX = Y and solve it
            startTime = tic;
            bd = reshape(blockedDofs{1}, [], 1);
            nbrDofs = size(X, 1);
            if isempty(bd)
                freeDofs = 1:nbrDofs;
            else
                freeDofs = ~any(1:nbrDofs == bd, 1);
            end
            partY = Y(freeDofs, :) - ...
                A(freeDofs, ~freeDofs)* ...
                X(~freeDofs, :);
            partA = A(freeDofs, freeDofs);
            partX = partA \ partY;
            X(freeDofs, :) = partX;
            fprintf("Solved system:\t\t\t\t%f secs\n", toc(startTime));
        end
    end
end

