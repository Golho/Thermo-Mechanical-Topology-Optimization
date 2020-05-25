classdef HeatFEM < HeatFEMBase
    %HEATFEM Summary of this class goes here
    %   Detailed explanation goes here
    
    properties
        spatialDimensions
    end
    
    properties(Access = protected)
        mesh % Gsmh object
        EnodElements % Mapping from Enod numbering to Gmsh numbering
    end
    
    methods
        function obj = HeatFEM(gmshData, elementType, tFinal, timeSteps, theta)
            %HEATFEM Construct an instance of this class
            %   Detailed explanation goes here
            if nargin < 5
                theta = 1;
            end
            obj.theta = theta;
            obj.mesh = gmshData;
            obj.spatialDimensions = elementType.dimensions;
            obj.nbrNodes = obj.mesh.nodes.numNodes;
            obj.nodeCoordinates = getGlobalCoordinates(obj.mesh);
            
            obj.ElementType = elementType;
            Enod = getElements(obj.mesh, obj.ElementType);
            obj.Enod = Enod(2:end, :);
            obj.EnodElements = Enod(1, :);
            x = obj.nodeCoordinates(1, :);
            y = obj.nodeCoordinates(2, :);
            z = obj.nodeCoordinates(3, :);
            obj.Ex = x(obj.Enod);
            obj.Ey = y(obj.Enod);
            obj.Ez = z(obj.Enod);
            
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
            filenameT = filePrefix + ".vtk";
            UnstructuredMat2VTK(filenameT, obj.nodeCoordinates, ...
                obj.Enod, obj.ElementType, fieldVector, "ascii", label, ...
                "NodeField");
        end
        
        function saveElementField(obj, filePrefix, fieldVector, label)
            filenameT = filePrefix + ".vtk";
            UnstructuredMat2VTK(filenameT, obj.nodeCoordinates, ...
                obj.Enod, obj.ElementType, fieldVector, "ascii", label, ...
                "ElementField", true);
        end
    end
    
    methods(Access = protected)
        function applyBoundaryConditions(obj)
            for i = 1:numel(obj.boundaryConditions)
                bc = obj.boundaryConditions{i};
                elementBlocks = getPhysicalEntity(obj.mesh, bc.physicalName);
                switch bc.type
                    case 'Neumann'
                        magnitude = bc.value;
                        func = @(coords, elementType, tag) obj.elementLoad(...
                            coords, elementType, magnitude);
                        startTime = tic;
                        f = obj.integrate(elementBlocks, func, 1);
                        fprintf("Integrated Neumann BC:\t\t%f secs\n", toc(startTime));
                        obj.fl(:, bc.timeSteps) = ...
                            obj.fl(:, bc.timeSteps) + ...
                            repmat(f, [1 length(bc.timeSteps)]);
                    case 'Robin'
                        magnitude = bc.alpha * bc.Tinf;
                        func = @(coords, elementType, tag) obj.elementLoad(...
                            coords, elementType, magnitude);
                        f = obj.integrate(elementBlocks, func, 1);
                        obj.fc(:, bc.timeSteps) = ...
                            obj.fc(:, bc.timeSteps) + ...
                            repmat(f, [1 length(bc.timeSteps)]);
                        
                        magnitude = bc.alpha;
                        func = @(coords, elementType, tag) obj.elementMass(...
                            coords, elementType, magnitude);
                        % TODO: handle different Kc at different time steps
                        obj.Kc = obj.Kc + obj.integrate(elementBlocks, func, 2);
                        
                    case 'Dirichlet'
                        startTime = tic;
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
                        fprintf("Applied Dirichlet BC:\t\t%f secs\n", toc(startTime));
                end
            end
        end
        
        function applyBodyConditions(obj)
            for i = 1:numel(obj.bodyConditions)
                bc = obj.bodyConditions{i};
                elementBlocks = getPhysicalEntity(obj.mesh, bc.physicalName);
                switch bc.type
                    case 'main'
                        startTime = tic;
                        func = @(coords, elementType, tag) obj.elementStiffness(...
                            coords, elementType, obj.material.Kappa);
                        obj.K = obj.K + obj.integrate(elementBlocks, func, 2);
                        
                        magnitude = obj.material.density * obj.material.heatCapacity;
                        func = @(coord, elementType, tag) obj.elementMass(...
                            coord, elementType, magnitude);
                        obj.M = obj.M + obj.integrate(elementBlocks, func, 2);
                        fprintf("Integrated main BC:\t\t\t%f secs\n", toc(startTime));
                    case 'load'
                        magnitude = bc.value;
                        func = @(coords, elementType, tag) obj.elementLoad(...
                            coords, elementType, magnitude);
                        f = obj.integrate(elementBlocks, func, 1);
                        obj.fv(:, bc.timeSteps) = ...
                            obj.fv(:, bc.timeSteps) + ...
                            repmat(f, [1 length(bc.timeSteps)]);
                    case 'initial'
                        startTime = tic;
                        % Extract all the nodes of the element blocks
                        nodes = [];
                        for elementBlock = elementBlocks
                            nodes = horzcat(nodes, [elementBlock.elements.nodeTags]);
                        end
                        bd = obj.Dofs(unique(nodes));
                        obj.temperatures(bd, 1) = bc.value;
                        fprintf("Applied initial BC:\t\t%f secs\n", toc(startTime));
                end
            end
        end
        
        function [matrix] = integrate(obj, elementBlocks, func, dim)
            numMatrixElements = 0;
            for block = elementBlocks
                numNodes = block.elementType.numNodes^dim * block.numElementsInBlock;
                numMatrixElements = numMatrixElements + numNodes;
            end
            I = zeros(numMatrixElements, 1);
            V = zeros(numMatrixElements, 1);
            if dim == 1
                J = 1;
            elseif dim == 2
                J = zeros(numMatrixElements, 1);
            else
                error('The input "dim" must be either 1 or 2');
            end
            counter = 0;
            for block = elementBlocks
                for element = block.elements
                    coord = obj.nodeCoordinates(:, element.nodeTags)';
                    edof = obj.Dofs(:, element.nodeTags);
                    elemMatrix = func(coord, block.elementType, element.elementTag);
                    sumDofs = numel(elemMatrix);
                    if dim == 1
                        I(counter+(1:sumDofs)) = edof;
                        V(counter+(1:sumDofs)) = elemMatrix;
                    elseif dim == 2
                        EdofEdof = repmat(edof(:), numel(edof), 1);
                        I(counter+(1:sumDofs)) = EdofEdof(:);
                        EEddooff = repelem(edof(:), numel(edof), 1);
                        J(counter+(1:sumDofs)) = EEddooff;
                        V(counter+(1:sumDofs)) = elemMatrix(:);
                    end
                    counter = counter + sumDofs;
                end
            end
            if dim == 1
                matrix = sparse(I(1:counter), J, V(1:counter), obj.nbrDofs, 1);
            elseif dim == 2
                matrix = sparse(I(1:counter), J(1:counter), V(1:counter), obj.nbrDofs, obj.nbrDofs);
            end
        end
    end
end
