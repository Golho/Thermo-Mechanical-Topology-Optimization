classdef MechFEMStructured < MechFEMBase
    %UNTITLED Summary of this class goes here
    %   Detailed explanation goes here
    
    properties
        spatialDimensions
    end
    
    properties(Access = protected)
        mesh
    end
    
    methods
        function obj = MechFEMStructured(structuredMesh, timeSteps, planarType)
            %HEATFEM Construct an instance of this class
            %   Detailed explanation goes here
            if nargin >= 3
                obj.planarType = planarType;
            end
            obj.timeSteps = timeSteps;
            
            obj.mesh = structuredMesh;
            obj.spatialDimensions = obj.mesh.Dimensions;
            obj.fieldDim = obj.spatialDimensions;
            obj.nbrNodes = obj.mesh.NumNodes;
            obj.nodeCoordinates = obj.mesh.coordinates();

            obj.init();
            
            obj.ElementType = obj.mesh.ElementType;
            [obj.Ex, obj.Ey, obj.Ez] = obj.mesh.elementCoordinates();
            obj.Enod = obj.mesh.elementNodes(1:obj.mesh.NumElements);
            obj.Edof = zeros(obj.ElementType.numNodes*obj.fieldDim, obj.mesh.NumElements);
            for dim = 1:obj.fieldDim
                dofs = obj.Dofs(dim, :);
                obj.Edof(dim:obj.fieldDim:end, :) = dofs(obj.Enod);
            end
        end
        
        function saveNodeField(obj, filePrefix, fieldVector, label)
            temp3DMatrix = zeros(obj.mesh.Nx, obj.mesh.Ny, obj.mesh.Nz);
            temp3DMatrix(:) = fieldVector;
            filenameT = filePrefix + ".vtk";
            Mat2VTK(filenameT, temp3DMatrix, "ascii", label, ...
                "NodeField");
        end
        
        function saveNodeVectorField(obj, filePrefix, fieldVector, label)
            % Make sure the padded field vector has zeros in the dimensions
            % not used
            paddedFieldVector = zeros(3, size(fieldVector, 2));
            paddedFieldVector(1:size(fieldVector, 1), :) = fieldVector;
            vector4DMatrix = zeros(3, obj.mesh.Nx, obj.mesh.Ny, obj.mesh.Nz);
            vector4DMatrix(:) = paddedFieldVector;
            filenameT = filePrefix + ".vtk";
            Vector2VTK(filenameT, vector4DMatrix, "ascii", label, ...
                "NodeField");
        end
        
        function saveElementField(obj, filePrefix, fieldVector, label)
            temp3DMatrix = zeros(max(obj.mesh.Nx-1, 1),...
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
        
        function f = getDummy(obj, name)
            f = sparse(obj.nbrDofs, obj.timeSteps);
            for iBc = 1:numel(obj.boundaryConditions)
                bc = obj.boundaryConditions{iBc};
                if bc.type == "dummy" && bc.name == name
                    assert(bc.value ~= 0, "The dummy load magnitude must be non-zero");
                    components = logical(reshape(bc.components, [], 1));
                    dofs = obj.Dofs(components, bc.nodes);
                    f(dofs, bc.timeSteps) = bc.value;
                end
            end
            if ~nnz(f)
                warning("The dummy load with the name %s was not found", name);
            end
        end
    end
    
    methods(Access = protected)
        function applyBoundaryConditions(obj)
            for i = 1:numel(obj.boundaryConditions)
                bc = obj.boundaryConditions{i};
                switch bc.type
                    case 'Neumann'
                        startTime = tic;
                        components = logical(reshape(bc.components, [], 1));
                        dofs = obj.Dofs(components, bc.nodes);
                        f = sparse(dofs, 1, bc.value, obj.nbrDofs, 1);
                        fprintf("Integrated Neumann BC:\t\t%f secs\n", toc(startTime));
                        obj.fl(:, bc.timeSteps) = ...
                            obj.fl(:, bc.timeSteps) + ...
                            repmat(f, [1 length(bc.timeSteps)]);
                    case 'Robin'
                        startTime = tic;
                        components = logical(reshape(bc.components, [], 1));
                        dofs = obj.Dofs(components, bc.nodes);
                        obj.K_c = sparse(dofs, dofs, bc.value, obj.nbrDofs, obj.nbrDofs);
                        fprintf("Applied Robin BC:\t\t%f secs\n", toc(startTime));
                    case 'Dirichlet'
                        startTime = tic;
                        % Extract all the nodes of the element blocks
                        components = logical(reshape(bc.components, [], 1));
                        nodes = bc.nodes;
                        bd = obj.Dofs(components, unique(nodes));
                        addBd = @(c) {[c; bd(:)]};
                        obj.blockedDofs(bc.timeSteps) = cellfun(addBd, ...
                            obj.blockedDofs(bc.timeSteps));
                        obj.displacements(bd(:), bc.timeSteps) = bc.value;
                        fprintf("Applied Dirichlet BC:\t\t%f secs\n", toc(startTime));
                    case "dummy"
                        continue
                    otherwise
                        warning("The boundary condition type %s is not supported", bc.type);
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
                        obj.K = obj.applyMain(elementCoord);
                        
                        thermalMatrix = obj.elementTempStiffness(elementCoord, ...
                            obj.ElementType, obj.material.stiffness, obj.material.thermalExp, obj.planarType);
                        obj.K_thermal = obj.integrate(thermalMatrix, 3);
                    case 'load'
                        loadMatrix = obj.elementLoad(elementCoord, obj.ElementType, bc.value);
                        f = obj.integrate(loadMatrix, 1);
                        obj.fv(:, bc.timeSteps) = ...
                            obj.fv(:, bc.timeSteps) + ...
                            repmat(f, [1 length(bc.timeSteps)]);
                end
            end
            obj.appliedBodies = true;
        end
        
        function [K] = applyMain(obj, elementCoord)
            startTime = tic;
            stiffMatrix = obj.elementStiffness(elementCoord, obj.ElementType, obj.material.stiffness, obj.planarType);
            K = obj.integrate(stiffMatrix, 2);
            fprintf("Integrated main BC:\t\t\t%f secs\n", toc(startTime));
        end
        
        function [matrix] = integrate(obj, elementMatrix, dim)
            nbrElements = obj.mesh.NumElements;
            V = repmat(elementMatrix(:), nbrElements, 1);
            if dim == 1
                % For vectors of size [nbrDofs x 1]
                I = obj.Edof(:);
                J = 1;
                matrix = sparse(I, J, V, obj.nbrDofs, 1);
            elseif dim == 2
                % For matrices of size [nbrDofs x nbrDofs]
                EdofEdof = repmat(obj.Edof, size(obj.Edof, 1), 1);
                I = EdofEdof(:);
                EEnnooff = repelem(obj.Edof(:), size(obj.Edof, 1), 1);
                J = EEnnooff;
                matrix = sparse(I, J, V, obj.nbrDofs, obj.nbrDofs);
            elseif dim == 3
                % For matrices of size [nbrDofs x nbrNodes], for example
                % the coupling between the displacements and temperatures
                EdofEdof = repmat(obj.Edof, size(obj.Enod, 1), 1);
                I = EdofEdof(:);
                EEnnoodd = repelem(obj.Enod(:), size(obj.Edof, 1), 1);
                J = EEnnoodd;
                matrix = sparse(I, J, V, obj.nbrDofs, obj.nbrNodes);
            end
        end
    end
end

