function [foundElementBlocks] = getPhysicalEntity(gmsh, physicalName)
% GETPHYSICALENTITY Extract the element entity blocks of all entity within
% a specified physical group


% Find the physicalName entity in the gmsh struct
found = 0;
for pn = gmsh.physicalNames.physicalNames
    if strcmp(pn.name, physicalName)
        disp('Found the entity');
        foundPn = pn;
        found = 1;
        break;
    end
end

if ~found
    error('No entity with the physical name "%s" given', physicalName);
end


% Decide which entity list to search for the physical tag based on the
% physical group dimension
switch foundPn.dimension
    case 0
        entities = gmsh.entities.points;
        tagField = 'pointTag';
    case 1
        entities = gmsh.entities.curves;
        tagField = 'curveTag';
    case 2
        entities = gmsh.entities.surfaces;
        tagField = 'surfaceTag';
    case 3
        entities = gmsh.entities.volumes;
        tagField = 'volumeTag';
    otherwise
        error('The dimension parameter of the physical group is invalid');
end

% Find in which entities the physical group is present
idx = arrayfun(@(entity) any(entity.physicalTags == foundPn.physicalTag), ...
    entities);
foundEntities = entities(idx);

compFun = @(entityBlock) any([foundEntities.(tagField)] == entityBlock.entityTag) && entityBlock.entityDim == foundPn.dimension;

idx = arrayfun(compFun, gmsh.elements.entityBlocks);
foundElementBlocks = gmsh.elements.entityBlocks(idx);

end

