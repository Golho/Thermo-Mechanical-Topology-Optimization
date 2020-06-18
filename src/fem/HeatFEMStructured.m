classdef HeatFEMStructured < HeatFEMBase
    %HEATFEM Summary of this class goes here
    %   Detailed explanation goes here
    
    properties
        spatialDimensions
        mesh
    end
    
    methods
        function obj = HeatFEMStructured(structuredMesh, tFinal, timeSteps, theta)
            %HEATFEM Construct an instance of this class
            %   Detailed explanation goes here
            if nargin < 4
                theta = 1;
            end
            obj.theta = theta;
            obj.mesh = structuredMesh;
            obj.spatialDimensions = obj.mesh.Dimensions;
            obj.nbrNodes = obj.mesh.NumNodes;
            obj.nodeCoordinates = obj.mesh.coordinates();
            
            [obj.Ex, obj.Ey, obj.Ez] = obj.mesh.elementCoordinates();
            obj.Enod = obj.mesh.elementNodes(1:obj.mesh.NumElements);
            obj.ElementType = obj.mesh.ElementType;
            
            if exist('tFinal', 'var')
                obj.transient = true;
                obj.tFinal = tFinal;
                if ~exist('timeSteps', 'var')
                    error(['The number of time steps must be specified if ' ...
                        'a transient model is chosem']);
                end
                obj.timeSteps = timeSteps;
            end
            
            obj.init();
        end
        
        function saveNodeField(obj, filePrefix, fieldVector, label)
            temp3DMatrix = zeros(...
                size(fieldVector, 1), ...
                obj.mesh.Nx, ...
                obj.mesh.Ny, ...
                obj.mesh.Nz);
            temp3DMatrix(:) = fieldVector;
            filenameT = filePrefix + ".vtk";
            Mat2VTK(filenameT, temp3DMatrix, "ascii", label, ...
                "NodeField");
        end
        
        function saveElementField(obj, filePrefix, fieldVector, label)
            temp3DMatrix = zeros(...
                size(fieldVector, 1), ...
                max(obj.mesh.Nx-1, 1), ...
                max(obj.mesh.Ny-1, 1), ...
                max(obj.mesh.Nz-1, 1));
            temp3DMatrix(:) = fieldVector;
            filenameT = filePrefix + ".vtk";
            filenamePointsT = filePrefix + "_POINTS.vtk";
            Mat2VTK(filenameT, temp3DMatrix, "ascii", label, ...
                "ElementField", true);
            Mat2VTK(filenamePointsT, temp3DMatrix, "ascii", label, ...
                "ElementField", false);
        end
    end
    
    methods(Access = protected)
        function applyBoundaryConditions(obj)
            for i = 1:numel(obj.boundaryConditions)
                bc = obj.boundaryConditions{i};
                switch bc.type
                    case 'Neumann'
                        startTime = tic;
                        f = sparse(bc.nodes, 1, bc.value, obj.nbrDofs, 1);
                        fprintf("Integrated Neumann BC:\t\t%f secs\n", toc(startTime));
                        obj.fl(:, bc.timeSteps) = ...
                            obj.fl(:, bc.timeSteps) + ...
                            repmat(f, [1 length(bc.timeSteps)]);
                    case 'Robin'
                        error("Not done implemented");
                        
                    case 'Dirichlet'
                        startTime = tic;
                        % Extract all the nodes of the element blocks
                        nodes = bc.nodes;
                        bd = obj.Dofs(unique(nodes));
                        addBd = @(c) {[c, bd]};
                        obj.blockedDofs(bc.timeSteps) = cellfun(addBd, ...
                            obj.blockedDofs(bc.timeSteps));
                        obj.temperatures(bd, bc.timeSteps) = bc.value;
                        fprintf("Applied Dirichlet BC:\t\t%f secs\n", toc(startTime));
                end
            end
            obj.appliedBoundaries = true;
        end
        
        function applyBodyConditions(obj)
            enod = obj.Enod(:, 1);
            elementCoord = obj.nodeCoordinates(:, enod)';
            for i = 1:numel(obj.bodyConditions)
                bc = obj.bodyConditions{i};
                switch bc.type
                    case 'main'
                        [obj.K, obj.M] = obj.applyMain(elementCoord);
                    case 'load'
                        loadMatrix = obj.elementLoad(elementCoord, bc.value, obj.thickness);
                        f = obj.integrate(loadMatrix, 1);
                        obj.fv(:, bc.timeSteps) = ...
                            obj.fv(:, bc.timeSteps) + ...
                            repmat(f, [1 length(bc.timeSteps)]);
                    case 'initial'
                        startTime = tic;
                        nodes = bc.nodes;
                        bd = obj.Dofs(unique(nodes));
                        obj.temperatures(bd, 1) = bc.value;
                        fprintf("Applied initial BC:\t\t%f secs\n", toc(startTime));
                end
            end
            obj.appliedBodies = true;
        end
        
        function [K, M] = applyMain(obj, elementCoord)
            startTime = tic;
            stiffMatrix = obj.elementStiffness(elementCoord, ...
                obj.mesh.ElementType, obj.material.Kappa);
            K = obj.integrate(stiffMatrix, 2);
            
            magnitude = obj.material.density * obj.material.heatCapacity;
            massMatrix = obj.elementMass(elementCoord, ...
                obj.mesh.ElementType, magnitude);
            M = obj.integrate(massMatrix, 2);
            fprintf("Integrated main BC:\t\t\t%f secs\n", toc(startTime));
        end
        
        function [matrix] = integrate(obj, elementMatrix, dim)
            nbrElements = obj.mesh.NumElements;
            V = repmat(elementMatrix(:), nbrElements, 1);
            if dim == 1
                I = obj.Enod(:);
                J = 1;
                matrix = sparse(I, J, V, obj.nbrDofs, 1);
            elseif dim == 2
                EnodEnod = repmat(obj.Enod, size(obj.Enod, 1), 1);
                I = EnodEnod(:);
                EEnnoodd = repelem(obj.Enod(:), size(obj.Enod, 1), 1);
                J = EEnnoodd;
                matrix = sparse(I, J, V, obj.nbrDofs, obj.nbrDofs);
            end
        end
    end
end
