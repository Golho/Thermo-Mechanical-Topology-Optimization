classdef heatFEM
    %HEATFEM Summary of this class goes here
    %   Detailed explanation goes here
    
    properties
        spatialDimensions
        nodeCoordinates
        edges
        elementData
        elementType
        % All the information about the mesh is given above
        % The properties below are just helpers
        nbrNodes
        nbrDofs
        Edof
        Enod
        Dof
        Ex
        Ey
        
        thickness
        D % conductivity matrix
        K % stiffness matrix
        blockedDofs
        
        temperatures
        loads
    end
    
    methods
        function obj = heatFEM(p, e, t, elementType)
            %HEATFEM Construct an instance of this class
            %   Detailed explanation goes here
            obj.nodeCoordinates = p';
            obj.edges = e;
            obj.elementData = t;
            obj.elementType = elementType;
            
            obj.nbrNodes = length(p);
            obj.nbrDofs = obj.nbrNodes;
            
            obj.Enod = zeros(size(t'));
            obj.Enod(:, 1) = 1:length(obj.Enod);
            obj.Enod(:, 2:end) = t(1:end-1, :)';
            obj.Dof = zeros(obj.nbrNodes, 1); % 1 dof per node for scalar field
            obj.Dof = 1:length(obj.Dof);
            
            [obj.Ex, obj.Ey] = coordxtr(obj.Enod, obj.nodeCoordinates);
            
            % As dof == node
            obj.Edof = obj.Enod;
            
            % Initialize result vectors
            obj.loads = zeros(obj.nbrDofs, 1);
            obj.temperatures = zeros(obj.nbrDofs, 1);
        end
        
        function obj = applyBoundaryConditions(obj, bc)
            obj.blockedDofs = bc(:, 1);
            obj.temperatures(bc(:, 1)) = bc(:, 2);
        end
        
        function obj = applyLoads(obj, f)
            obj.loads = f;
        end
        
        function obj = setMaterial(obj, material, thickness)
            obj.D = material.D;
            obj.thickness = thickness;
            % Add other properties like density and heat capacity for
            % transient problem
        end
        
        function obj = assembleStiffness(obj)
            %METHOD1 Summary of this method goes here
            %   Detailed explanation goes here
            for e = 1:length(obj.Edof)
                edof = obj.Edof(e, 2:end);
                ex = obj.Ex(e, 2:end);
                ey = obj.Ey(e, 2:end);
                switch obj.elementType
                    case 3 % QUADS
                        Ke = flw2re(ex, ey, obj.D, obj.thickness);
                    otherwise
                        error("The element stiffness matrix is not yet implemented for the current element type");
                end
                obj.K(edof, edof) = obj.K(edof, edof) + Ke;
            end
        end
        
        function obj = solve(obj)
            freeDofs = setdiff(obj.blockedDofs, 1:obj.nbrDofs);
            loadsPartitioned = obj.loads(freeDofs) - ...
                obj.K(freeDofs, obj.blockedDofs);
            KPartitioned = obj.K(freeDofs, freeDofs);
            
            solution = KPartitioned\loadsPartitioned;
            obj.temperatures(freeDofs) = solution;
        end
    end
end

