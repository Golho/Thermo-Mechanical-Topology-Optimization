function [Enod] = getElements(gmshData, elementType)
%GETGLOBALNODES Summary of this function goes here
%   Detailed explanation goes here
Enod = [];
for entityBlock = gmshData.elements.entityBlocks
    if entityBlock.elementType == elementType
        elements = [entityBlock.elements.elementTag];
        nodes = reshape([entityBlock.elements.nodeTags], [], length(elements)); 
        enod = [elements; nodes];
        Enod = horzcat(Enod, enod);
    end
end
end

