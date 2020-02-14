function msh = gmshParser(filename)

%% Reads a mesh in msh format, version 4.1
%
% Usage: 
% To define all variables m.LINES, M.TRIANGLES, etc 
% (Warning: This creates a very large structure. Do you really need it?)
%            m = load_gmsh4('a.msh')
%
% To define only certain variables (for example TETS and HEXS)
%            m = load_gmsh4('a.msh', [ 5 6])
%
% To define no variables (i.e. if you prefer to use m.ELE_INFOS(i,2))
%            m = load_gmsh4('a.msh', -1)
%            m = load_gmsh4('a.msh', [])
%
% Copyright (C) 2007  JP Moitinho de Almeida (moitinho@civil.ist.utl.pt)
% and  R Lorphevre(r(point)lorphevre(at)ulg(point)ac(point)be)
%
% based on load_gmsh.m supplied with gmsh-2.0
%
% Structure msh always has the following elements:
%
% msh.MIN, msh.MAX - Bounding box
% msh.nbNod - Number of nodes
% msh.nbElm - Total number of elements
% msh.nbType(i) - Number of elements of type i (i in 1:19)
% msh.POS(i,j) - j'th coordinate of node i (j in 1:3, i in 1:msh.nbNod)
% msh.ELE_INFOS(i,1) - id of element i (i in 1:msh.nbElm)
% msh.ELE_INFOS(i,2) - type of element i
% msh.ELE_INFOS(i,3) - number of tags for element i
% msh.ELE_INFOS(i,4) - Dimension (0D, 1D, 2D or 3D) of element i
% msh.ELE_TAGS(i,j) - Tags of element i (j in 1:msh.ELE_INFOS(i,3))
% msh.ELE_NODES(i,j) - Nodes of element i (j in 1:k, with
%                       k = msh.NODES_PER_TYPE_OF_ELEMENT(msh.ELE_INFOS(i,2)))
%
% These elements are created if requested:
%
% msh.nbLines - number of 2 node lines
% msh.LINES(i,1:2) - nodes of line i
% msh.LINES(i,3) - tag (WHICH ?????) of line i
%
% msh.nbTriangles - number of 2 node triangles
% msh.TRIANGLES(i,1:3) - nodes of triangle i
% msh.TRIANGLES(i,4) - tag (WHICH ?????) of triangle i
%
% Etc

% These definitions need to be updated when new element types are added to gmsh
%
% msh.Types{i}{1} Number of an element of type i
% msh.Types{i}{2} Dimension (0D, 1D, 2D or 3D) of element of type i
% msh.Types{i}{3} Name to add to the structure with elements of type i
% msh.Types{i}{4} Name to add to the structure with the number of elements of type i
%

msh.elementTypes = { ...
    { 2, 1, 'LIN_2', 'nbLines'}, ... % 1
    { 3,  2, 'TRI_3', 'nbTriangles'}, ...
    { 4,  2, 'QUA_4', 'nbQuads'}, ...  
    { 4,  3, 'TET_4', 'nbTets'}, ...
    { 8,  3, 'HEX_8', 'nbHexas'}, ... %5
    { 6,  3, 'PRI_6', 'nbPrisms'}, ...
    { 5,  3, 'PYR_5', 'nbPyramids'}, ...
    { 3,  1, 'LIN_3', 'nbLines3'}, ...
    { 6,  2, 'TRI_6', 'nbTriangles6'}, ...
    { 9,  2, 'QUA_9', 'nbQuads9'}, ... % 10
    { 10,  3, 'TET_10', 'nbTets10'}, ...
    { 27,  3, 'HEX_27', 'nbHexas27'}, ...
    { 18,  3, 'PRI_18', 'nbPrisms18'}, ...
    { 14,  3, 'PYR_14', 'nbPyramids14'}, ...
    { 1,  0, 'PNT', 'nbPoints'}, ... % 15
    { 8,  3, 'QUA_8', 'nbQuads8'}, ...
    { 20,  3, 'HEX_20', 'nbHexas20'}, ...
    { 15,  3, 'PRI_15', 'nbPrisms15'}, ...
    { 13,  3, 'PYR_13', 'nbPyramids13'}, ...
};
% TODO: Add the new element types of GMSH, see 
% https://gitlab.onelab.info/gmsh/gmsh/blob/master/Common/GmshDefines.h
% TODO: Add error handling for input parsing in sub functions

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
            msh.elements = readElements(fid, msh.elementTypes);
    end
    tline = fgetl(fid);
end

fclose(fid);

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

function [elements] = readElements(fid, elementTypes)
    metaData = fscanf(fid, '%f', 4); % <numEntityBlocks> <numElements> <minElementTag> <maxElementTag>
    
    elements.numEntityBlocks = metaData(1);
    elements.numElements= metaData(2);
    elements.minElementTag = metaData(3);
    elements.maxElementTag = metaData(4);
    
    for iEntity = 1:elements.numEntityBlocks
        entityData = fscanf(fid, '%d', 4); % <entityDim> <entityTag> <elementType> <numNodesInBlock>
        
        numNodes = elementTypes{entityData(3)}{1};
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
            'elementType', entityData(3), ...
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