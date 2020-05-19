classdef Elements
    %ELEMENTS Summary of this class goes here
    %   Detailed explanation goes here
    properties
        id
        numNodes
        dimensions
        vtkCellType
    end
    
    methods
        function e = Elements(id, numNodes, dimensions, vtkCellType)
            e.id = id;
            e.numNodes = numNodes;
            e.dimensions = dimensions;
            if nargin < 4
                e.vtkCellType = 0; % Value of 0 indicate non-existense in VTK
            else
                e.vtkCellType = vtkCellType;
            end
        end
    end
    
    methods(Static)
        function [enum] = getEnum(id)
            enums = enumeration("Elements");
            for enum = enums'
                if enum.id == id
                    return;
                end
            end
            error("Could not find enum");
        end
    end
    
    enumeration
        LIN_2   (1, 2, 1, 3)
        TRI_3   (2, 3, 2, 5)
        QUA_4   (3, 4, 2, 9)
        TET_4   (4, 4, 3, 10)
        HEX_8   (5, 8, 3, 12)
        PRI_6   (6, 6, 3, 13)
        PYR_5   (7, 5, 3, 14)
        LIN_3   (8, 3, 1)
        TRI_6   (9, 6, 2)
        QUA_9   (10, 9, 2)
        TET_10  (11, 10, 3)
        HEX_27  (12, 27, 3)
        PRI_18  (13, 18, 3)
        PYR_14  (14, 14, 3)
        PNT     (15, 1, 0, 1)
    end 
end

