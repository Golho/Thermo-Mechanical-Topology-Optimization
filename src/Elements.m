classdef Elements
    %ELEMENTS Summary of this class goes here
    %   Detailed explanation goes here
    properties
        id
        numNodes
        dimensions
    end
    
    methods
        function e = Elements(id, numNodes, dimensions)
            e.id = id;
            e.numNodes = numNodes;
            e.dimensions = dimensions;
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
        LIN_2   (1, 2, 1)
        TRI_3   (2, 3, 2)
        QUA_4   (3, 4, 2)
        TET_4   (4, 4, 3)
        HEX_8   (5, 8, 3)
        PRI_6   (6, 6, 3)
        PYR_5   (7, 5, 3)
        LIN_3   (8, 3, 1)
        TRI_6   (9, 6, 2)
        QUA_9   (10, 9, 2)
        TET_10  (11, 10, 3)
        HEX_27  (12, 27, 3)
        PRI_18  (13, 18, 3)
        PYR_14  (14, 14, 3)
        PNT     (15, 1, 0)
    end 
end

