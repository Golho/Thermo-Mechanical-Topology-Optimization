classdef FEMBase < matlab.mixin.Copyable
    %UNTITLED Summary of this class goes here
    %   Detailed explanation goes here
    
    properties(Abstract)
        configuration
        fieldDim % Dimensions of the field variable (> 1 = vector field)
        spatialDimensions
    end
    
    properties
        nbrNodes
        nbrDofs
        timeSteps = 0;
        tFinal = -1;
        
        % Main elements properties
        Enod
        Edof
        Ex
        Ey
        Ez
        ElementType
    end
    
    properties(Access = protected)
        Dofs
        transient = false;
        material = Material();
        boundaryConditions
        bodyConditions
        
        appliedBoundaries = false;
        appliedBodies = false;
    end
    
    methods(Abstract)
        assemble(obj);
        solve(obj);
        saveNodeField(obj, filePrefix, fieldVector);
        saveElementField(obj, filePrefix, fieldVector);
    end
    
    methods
        function setMaterial(obj, material)
            obj.material = material;
        end
        
        function addBoundaryCondition(obj, boundaryCondition)
            obj.boundaryConditions{end + 1} = boundaryCondition;
        end
        
        function addBodyCondition(obj, bodyCondition)
            obj.bodyConditions{end + 1} = bodyCondition;
        end
        
        function nbrDofs = get.nbrDofs(obj)
            nbrDofs = obj.nbrNodes * obj.fieldDim;
        end
        
        function dofs = getDofs(obj, nodes)
            dofs = obj.Dofs(:, nodes);
            dofs = dofs(:);
        end
    end
end

