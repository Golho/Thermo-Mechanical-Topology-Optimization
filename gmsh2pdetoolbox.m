%GMSH2PDETOOLBOX Reads a mesh in msh format, version 4.1 and returns the
%arrays p, e and t according to the MATLAB pde toolbox syntax. Filename
%refers to the .msh file. Note that only 2D meshes are processed
%since these are the only elements allowed in pde toolbox.
%
%SEE ALSO load_gmsh load_gmsh4 initmesh
%
% Copyright (C) 2014 Arthur Levy. https://github.com/arthurlevy/
%
% This fucntion uses the routine load_gmsh4 in load_gmsh2.m: Copyright (C)
% 2007  JP Moitinho de Almeida (moitinho@civil.ist.utl.pt) and  R
% Lorphevre(r(point)lorphevre(at)ulg(point)ac(point)be) It is based on
% load_gmsh.m supplied with gmsh-2.0
function [p, e, t] = gmsh2pdetoolbox(filename, spatialDimension, elementType)

%loads the gmsh file using JP Moitinho de Almeida (moitinho@civil.ist.utl.pt)
% and  R Lorphevre(r(point)lorphevre(at)ulg(point)ac(point)be) code.
meshStructure = loadGmsh4(filename,-1);


%% NODES
%nodes positions
disp('importing nodes');
% p should have dimensions [spatialDimension x numNodes]
p = meshStructure.POS(:,1:spatialDimension)';


%% EDGES
% TODO: Add support for edges
% disp('importing edges');
% %find the edges (ie dimension 1)
% is_edge = (meshStructure.ELE_INFOS(:,2)==1);
% %add edge data
% e = meshStructure.ELE_NODES(is_edge,1:2)';
% e(3,:) = 0;
% e(4,:) = 1;
% %tag is important for applying boundary conditions
% e(5,:) = meshStructure.ELE_TAGS(is_edge,1)';
% e(6,:) = 1;
% e(7,:) = 0;
e = [];

%% TRIANGLES
sprintf('importing elements: %s ', meshStructure.Types{elementType}{3})
% count the number of elements with the specified element type
numElements = 0;
for entity = meshStructure.ELE_DATA'
    if entity.elementType == elementType
        numElements = numElements + length(entity.elements);
    end
end
numNodesPerElements = meshStructure.Types{elementType}{1};
t = zeros(numNodesPerElements + 1, numElements);
% Throw away the element numbering from GMSH, as the output for Matlab is
% not dependent on it (and it doesn't have to start from 1 to numElements)
i = 0;
for entity = meshStructure.ELE_DATA'
    if entity.elementType == elementType
        elementsInEntity = length(entity.elements);
        t(1:end-1, (1:elementsInEntity) + i) = entity.elements(:, 2:end)';
    end
end
t(end,:) = 1;

end


function msh = loadGmsh4(filename, which)

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

nargchk(1, 2, nargin);

msh.Types = { ...
    { 2, 1, 'LINES', 'nbLines'}, ... % 1
    { 3,  2, 'TRIANGLES', 'nbTriangles'}, ...
    { 4,  2, 'QUADS', 'nbQuads'}, ...  
    { 4,  3, 'TETS', 'nbTets'}, ...
    { 8,  3, 'HEXAS', 'nbHexas'}, ... %5
    { 6,  3, 'PRISMS', 'nbPrisms'}, ...
    { 5,  3, 'PYRAMIDS', 'nbPyramids'}, ...
    { 3,  1, 'LINES3', 'nbLines3'}, ...
    { 6,  2, 'TRIANGLES6', 'nbTriangles6'}, ...
    { 9,  2, 'QUADS9', 'nbQuads9'}, ... % 10
    { 10,  3, 'TETS10', 'nbTets10'}, ...
    { 27,  3, 'HEXAS27', 'nbHexas27'}, ...
    { 18,  3, 'PRISMS18', 'nbPrisms18'}, ...
    { 14,  3, 'PYRAMIDS14', 'nbPyramids14'}, ...
    { 1,  0, 'POINTS', 'nbPoints'}, ... % 15
    { 8,  3, 'QUADS8', 'nbQuads8'}, ...
    { 20,  3, 'HEXAS20', 'nbHexas20'}, ...
    { 15,  3, 'PRISMS15', 'nbPrisms15'}, ...
    { 13,  3, 'PYRAMIDS13', 'nbPyramids13'}, ...
};
% TODO: Add the new element types of GMSH, see 
% https://gitlab.onelab.info/gmsh/gmsh/blob/master/Common/GmshDefines.h
% TODO: Add error handling for input parsing in sub functions
                     
ntypes = length(msh.Types);

if (nargin==1)
    which = 1:ntypes;
else
    if isscalar(which) && which == -1
        which = [];
    end
end

% Could check that "which" is properly defined....

fid = fopen(filename, 'r');
fileformat = 0; % Assume wrong file

tline = fgetl(fid);
if (feof(fid))
    fprintf('Empty file: %s \n',  filename);
    msh = [];
    return;
end

%% Read mesh type
if (strcmp(tline, '$MeshFormat'))
    fileformat = 2;
    tline = fgetl(fid);
    if (feof(fid))
        fprintf('Syntax error (no version) in: %s\n',  filename);
        fileformat = 0;
    end
    [ form ] = sscanf( tline, '%f %d %d'); % <version> <file-type> <data-size>
    if (form(1) ~= 4.1)
        fprintf('The GMSH file version %f, is not supported (only version 4.1 is) \n', form(1) );
        fileformat = 0;
    end
    if (form(2) ~= 0)
        disp('Sorry, this program can only read ASCII format');
        fileformat = 0;
    end
    fgetl(fid);    % this should be $EndMeshFormat
    if (feof(fid))
        fprintf('Syntax error (no $EndMeshFormat) in: %s \n',  filename);
        fileformat = 0;
    end
end

if (~fileformat)
    msh = [];
    return
end

%% Read the remaining lines and find the blocks with $Nodes and $Elements
while ~feof(fid)
    tline = fgetl(fid);
    switch tline
        case '$Nodes'
            msh.POS = readNodeData(fid);
            msh.nbNod = length(msh.POS);
            msh.MAX = max(msh.POS);
            msh.MIN = min(msh.POS);
        case '$Elements'
            [msh.ELE_HEADER, msh.ELE_DATA] = readElementData(fid, msh.Types);
    end
end

if (fileformat == 0)
    msh = [];
end

fclose(fid);

end

function [positions] = readNodeData(fid)
    metaData = fscanf(fid, '%f', 4); % <numEntityBlocks> <numNodes> <minNodeTag> <maxNodeTag>
    positions = zeros(metaData(4), 3);
    fileformat = 1;
    for iEntityBlock = 1:metaData(1)
        entityData = fscanf(fid, '%d', 4); % <entityDim> <entityTag> <parametric> <numNodesInBlock>
        if entityData(3)
            disp('The parametric option is not yet supported by the parser');
            fileformat = 0;
            break;
        end
        tags = fscanf(fid, '%g', entityData(4));
        pos = fscanf(fid, '%g', [3, entityData(4)]); % transpose pos as fscanf inserts in column order
        positions(tags, :) = pos';
    end
    fgetl(fid);     % flush pointer to new line
    tline = fgetl(fid);     % this should be $EndNodes
    if (~strcmp(tline, '$EndNodes'))
        fprintf('Syntax error (no $EndNodes) in: %s \n',  filename);
        fileformat = 0;
    end
end

function [header, entities] = readElementData(fid, elementTypes)
    metaData = fscanf(fid, '%f', 4); % <numEntityBlocks> <numElements> <minElementTag> <maxElementTag>
    header = metaData;
    entity.entityDim = 0;
    entity.entityTag = 0;
    entity.elementType = 0;
    entity.numElementsInBlock = 0;
    entity.elements = [];
    entities = repmat(entity, [metaData(1), 1]);
    for iEntity = 1:metaData(1)
        entityData = fscanf(fid, '%d', 4); % <entityDim> <entityTag> <elementType> <numNodesInBlock>
        entities(iEntity).entityDim = entityData(1);
        entities(iEntity).entityTag = entityData(2);
        entities(iEntity).elementType = entityData(3);
        entities(iEntity).numElementsInBlock = entityData(4);
        numNodes = elementTypes{entityData(3)}{1};
         % transpose elements as fscanf inserts in column order
        entities(iEntity).elements = fscanf(fid, '%f', [1+numNodes entityData(4)])'; % <elementTag> <nodeTag> ... 
    end
    
    fgetl(fid);     % flush pointer to new line
    tline = fgetl(fid);     % this should be $EndElements
    if (~strcmp(tline, '$EndElements'))
        fprintf('Syntax error (no $EndElements) in: %s \n',  filename);
        fileformat = 0;
    end
end