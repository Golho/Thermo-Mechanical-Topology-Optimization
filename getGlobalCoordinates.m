function [nodeCoordinates] = getGlobalCoordinates(gmshData)
%GETGLOBALNODES Summary of this function goes here
%   Detailed explanation goes here
nodeCoordinates = zeros(gmshData.nodes.numNodes, 3);
for entityBlock = gmshData.nodes.entityBlocks
    nodeCoordinates(entityBlock.nodeTags, 1) = entityBlock.x;
    nodeCoordinates(entityBlock.nodeTags, 2) = entityBlock.y;
    nodeCoordinates(entityBlock.nodeTags, 3) = entityBlock.z;
end
end

