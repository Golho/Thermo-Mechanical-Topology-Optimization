function [nodeCoordinates] = getGlobalCoordinates(gmshData)
%GETGLOBALNODES Summary of this function goes here
%   Detailed explanation goes here
nodeCoordinates = zeros(3, gmshData.nodes.numNodes);
for entityBlock = gmshData.nodes.entityBlocks
    nodeCoordinates(1, entityBlock.nodeTags) = entityBlock.x;
    nodeCoordinates(2, entityBlock.nodeTags) = entityBlock.y;
    nodeCoordinates(3, entityBlock.nodeTags) = entityBlock.z;
end
end

