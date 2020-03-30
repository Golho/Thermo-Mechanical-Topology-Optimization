function msh = gmshParser(filename)

%% Reads a mesh in msh format, version 4.1

% TODO: Add the new element types of GMSH, see 
% https://gitlab.onelab.info/gmsh/gmsh/blob/master/Common/GmshDefines.h
% TODO: Add error handling for input parsing in sub functions

startTime = tic;

fid = fopen(filename, 'r');

tline = fgetl(fid);
if (feof(fid))
    fprintf('Empty file: %s \n',  filename);
    return;
end

%% Read the remaining lines and find the blocks with $Nodes and $Elements
while ~feof(fid)
    switch tline
        case '$MeshFormat' % Mandatory block
            msh.meshFormat = readMeshFormat(fid);
        case '$PhysicalNames'
            msh.physicalNames = readPhysicalNames(fid);
        case '$Entities'
            msh.entities = readEntities(fid);
        case '$Nodes'
            msh.nodes = readNodes(fid);
        case '$Elements'
            msh.elements = readElements(fid);
    end
    tline = fgetl(fid);
end

fclose(fid);
elapsed = toc(startTime);
fprintf("Parsed: %d nodes, %d elements in %f secs\n", msh.nodes.numNodes, ...
    msh.elements.numElements, elapsed)

end

function [meshFormat] = readMeshFormat(fid)
    data = fscanf(fid, '%f %d %d', 3);
    meshFormat.version = data(1);
    meshFormat.fileType = data(2);
    meshFormat.dataSize = data(3);
    
    if (meshFormat.version ~= 4.1)
        error('The GMSH file version %f, is not supported (only version 4.1 is) \n', form(1) );
    end

    if (meshFormat.fileType ~= 0)
        error('Sorry, this program can only read ASCII format');
    end
    
    flushBlock(fid, '$EndMeshFormat');
end

function [entities] = readEntities(fid)
    numData = fscanf(fid, '%d', 4);
    entities.numPoints = numData(1);
    entities.numCurves = numData(2);
    entities.numSurfaces = numData(3);
    entities.numVolumes = numData(4);
    
    for iPoint = 1:entities.numPoints
        data = fscanf(fid, '%d %f %f %f', 4);
        numPhysicalTags = fscanf(fid, '%d', 1);
        physicalTags = fscanf(fid, '%d', numPhysicalTags);
        entities.points(iPoint) = struct(...
            'pointTag', data(1), ...
            'X', data(2), ...
            'Y', data(3), ...
            'Z', data(4), ...
            'numPhysicalTags', numPhysicalTags, ...
            'physicalTags', physicalTags ...
        );
    end
    
    for iCurve = 1:entities.numCurves
        data = fscanf(fid, '%d %f %f %f %f %f %f', 7);
        numPhysicalTags = fscanf(fid, '%d', 1);
        physicalTags = fscanf(fid, '%d', numPhysicalTags);
        numBoundingPoints = fscanf(fid, '%d', 1);
        pointTags = fscanf(fid, '%d', numBoundingPoints);
        entities.curves(iCurve) = struct(...
            'curveTag', data(1), ...
            'minX', data(2), ...
            'minY', data(3), ...
            'minZ', data(4), ...
            'maxX', data(5), ...
            'maxY', data(6), ...
            'maxZ', data(7), ...
            'numPhysicalTags', numPhysicalTags, ...
            'physicalTags', physicalTags, ...
            'numBoundingPoints', numBoundingPoints, ...
            'pointTags', pointTags ...
        );
    end
    
    for iSurface = 1:entities.numSurfaces
        data = fscanf(fid, '%d %f %f %f %f %f %f', 7);
        numPhysicalTags = fscanf(fid, '%d', 1);
        physicalTags = fscanf(fid, '%d', numPhysicalTags);
        numBoundingCurves = fscanf(fid, '%d', 1);
        curveTags = fscanf(fid, '%d', numBoundingCurves);
        entities.surfaces(iSurface) = struct(...
            'surfaceTag', data(1), ...
            'minX', data(2), ...
            'minY', data(3), ...
            'minZ', data(4), ...
            'maxX', data(5), ...
            'maxY', data(6), ...
            'maxZ', data(7), ...
            'numPhysicalTags', numPhysicalTags, ...
            'physicalTags', physicalTags, ...
            'numBoundingCurves', numBoundingCurves, ...
            'curveTags', curveTags ...
        );
    end
    
    for iVolume = 1:entities.numVolumes
        data = fscanf(fid, '%d %f %f %f %f %f %f', 7);
        numPhysicalTags = fscanf(fid, '%d', 1);
        physicalTags = fscanf(fid, '%d', numPhysicalTags);
        numBoundingSurfaces = fscanf(fid, '%d', 1);
        surfaceTags = fscanf(fid, '%d', numBoundingSurfaces);
        entities.volumes(iVolume) = struct(...
            'volumeTag', data(1), ...
            'minX', data(2), ...
            'minY', data(3), ...
            'minZ', data(4), ...
            'maxX', data(5), ...
            'maxY', data(6), ...
            'maxZ', data(7), ...
            'numPhysicalTags', numPhysicalTags, ...
            'physicalTags', physicalTags, ...
            'numBoundingSurfaces', numBoundingSurfaces, ...
            'surfaceTags', surfaceTags ...
        );
    end
    
    flushBlock(fid, '$EndEntities');
end

function [pNames] = readPhysicalNames(fid)
    pNames.numPhysicalNames = fscanf(fid, '%d', 1);
    % Use textscan instead of fscanf to parse string
    data = textscan(fid, '%d %d %s', [3 pNames.numPhysicalNames])';
    % Fill structure array backwards to hopefully allocate all memory at
    % once
    for ips = pNames.numPhysicalNames:-1:1
        name = data{3}{ips};
        % trim the quotes from the start and the end of the string
        name = name(2:end-1);
        
        pNames.physicalNames(ips) = struct(...
            'dimension', data{1}(ips), ...
            'physicalTag', data{2}(ips), ...
            'name', name ...
        );
    end

    flushBlock(fid, '$EndPhysicalNames');
end

function [nodes] = readNodes(fid)
    metaData = fscanf(fid, '%d', 4); % <numEntityBlocks> <numNodes> <minNodeTag> <maxNodeTag>
    nodes.numEntityBlocks = metaData(1);
    nodes.numNodes = metaData(2);
    nodes.minNodeTag = metaData(3);
    nodes.maxNodeTag = metaData(4);
    
    for iEntityBlock = 1:nodes.numEntityBlocks
        entityData = fscanf(fid, '%d', 4); % <entityDim> <entityTag> <parametric> <numNodesInBlock>
        if entityData(3)
            disp('The parametric option is not yet supported by the parser');
            break;
        end
        nodeTags = fscanf(fid, '%g', entityData(4));
        pos = fscanf(fid, '%g', [3, entityData(4)])'; % transpose pos as fscanf inserts in column order
        u = []; v = []; w = [];  
        
        nodes.entityBlocks(iEntityBlock) = struct(...
            'entityDim', entityData(1), ...
            'entityTag', entityData(2), ...
            'parametric', entityData(3), ...
            'numNodesInBlock', entityData(4), ...
            'nodeTags', nodeTags, ...
            'x', pos(:, 1), ...
            'y', pos(:, 2), ...
            'z', pos(:, 3), ...
            'u', u, ...
            'v', v, ...
            'w', w ...
        );
    end
    
    flushBlock(fid, '$EndNodes');
end

function [elements] = readElements(fid)
    metaData = fscanf(fid, '%f', 4); % <numEntityBlocks> <numElements> <minElementTag> <maxElementTag>
    
    elements.numEntityBlocks = metaData(1);
    elements.numElements= metaData(2);
    elements.minElementTag = metaData(3);
    elements.maxElementTag = metaData(4);
    
    for iEntity = 1:elements.numEntityBlocks
        entityData = fscanf(fid, '%d', 4); % <entityDim> <entityTag> <elementType> <numNodesInBlock>
        elementType = Elements.getEnum(entityData(3));
        numNodes = elementType.numNodes;
         % transpose elements as fscanf inserts in column order
        elementDataMatrix = fscanf(fid, '%f', [1+numNodes entityData(4)])'; % <elementTag> <nodeTag> ... 
        
        for iElem = size(elementDataMatrix, 1):-1:1
            elems(iElem) = struct( ...
                'elementTag', elementDataMatrix(iElem, 1), ...
                'nodeTags', elementDataMatrix(iElem, 2:end) ...
            );
        end
        elements.entityBlocks(iEntity) = struct( ...
            'entityDim', entityData(1), ...
            'entityTag', entityData(2), ...
            'elementType', elementType, ...
            'numElementsInBlock', entityData(4), ...
            'elements', elems ...
        );
    end
    
    flushBlock(fid, '$EndElements');
end

function flushBlock(fid, endTag)
    fgetl(fid);     % flush pointer to new line
    tline = fgetl(fid);     % this should be $EndElements
    if (~strcmp(tline, endTag))
        fprintf('Syntax error (no %s) \n',  endTag);
    end
end