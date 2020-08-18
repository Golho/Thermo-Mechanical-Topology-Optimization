classdef VTKDumper < handle
    %UNTITLED Summary of this class goes here
    %   Detailed explanation goes here
    
    properties
        headerName
        meshData
        structured
        elementType
        pointData = struct(...
            "name", [], ...
            "type", [], ...
            "data", [], ...
            "nbrComponents", [], ...
            "timeSteps", [] ...
            );
        nbrPointData = 0;
        cellData = struct(...
            "name", [], ...
            "type", [], ...
            "data", [], ...
            "nbrComponents", [], ...
            "timeSteps", [] ...
            );
        nbrCellData = 0;
        maxTimeStep = 1;
        
        extrudeZ = false;
    end
    
    methods
        function obj = VTKDumper(headerName, meshData, elementType)
            if isa(meshData, "Gmsh")
                assert(nargin == 3, "The element type must be specified for a GMSH mesh");
                obj.structured = false;
                obj.elementType = elementType;
            elseif isa(meshData, "StructuredMesh")
                assert(nargin ~= 3, "The element type is supertfluous when mesh is structured");
                obj.structured = true;
            else
                error("The identified mesh type is not supported");
            end
            obj.headerName = headerName;
            obj.meshData = meshData;
        end
        
        function addData(obj, fieldData, dataType)
            %METHOD1 Summary of this method goes here
            %   Detailed explanation goes here
            % The number of entries in fieldData.data must equal the number
            % of cell/points times the number of components
            if dataType == "cell"
                if obj.structured
                    nbrDataEntities = obj.meshData.NumElements;
                else
                    nbrDataEntities = size(getElements(obj.meshData), 2);
                end
                obj.cellData(obj.nbrCellData+1) = fieldData;
                obj.nbrCellData = obj.nbrCellData + 1;
            elseif dataType == "point"
                if obj.structured
                    nbrDataEntities = obj.meshData.NumNodes;
                else
                    nbrDataEntities = obj.meshData.numNodes;
                end
                obj.pointData(obj.nbrPointData + 1) = fieldData;
                obj.nbrPointData = obj.nbrPointData + 1;
            end
            assert(numel(fieldData.data) == nbrDataEntities * fieldData.nbrComponents, ...
                "The length of the data is incompatible with the mesh");
            if fieldData.type == "vectors"
                assert(fieldData.nbrComponents == 3, "The vector field must have 3 components");
            end
            obj.maxTimeStep = max(max(fieldData.timeSteps), obj.maxTimeStep);
        end
        
        function dump(obj, filePrefix)
            for timeStep = 1:obj.maxTimeStep
                filename = filePrefix + "_" + timeStep + ".vtk";
                % Open the file.
                fid = fopen(filename, "w");
                if fid == -1
                    error("Cannot open file for writing.");
                end
                
                if obj.structured
                    if obj.meshData.Nz == 1
                        % Extend 2D geometries to be able to use volumes in
                        % Paraview
                        obj.meshData.Nz = 2;
                        obj.meshData.Lz = min(obj.meshData.Lx, obj.meshData.Ly);
                        obj.extrudeZ = true;
                        for iField = 1:obj.nbrPointData
                            obj.pointData(iField).data = ...
                                repmat(obj.pointData(iField).data, 1, 2);
                        end
                    end
                    obj.writeStructuredHeader(fid);
                else
                    obj.writeGmshHeader(fid);
                end
                
                if numel(obj.cellData) > 0
                    fprintf(fid, "CELL_DATA %d\n", obj.meshData.NumElements);
                    for fieldData = obj.cellData
                        if ismember(timeStep, fieldData.timeSteps)
                            obj.writeDataset(fid, fieldData);
                        end
                    end
                end
                
                if numel(obj.cellData) > 0
                    fprintf(fid, "POINT_DATA %d\n", obj.meshData.NumNodes);
                    for fieldData = obj.pointData
                        if ismember(timeStep, fieldData.timeSteps)
                            obj.writeDataset(fid, fieldData);
                        end
                    end
                end
                
                fclose(fid);
            end
        end
    end
    
    methods(Access = protected)
        function writeStructuredHeader(obj, fileId)
            dx = obj.meshData.Lx / max(obj.meshData.Nx - 1, 1);
            dy = obj.meshData.Ly / max(obj.meshData.Ny - 1, 1);
            dz = obj.meshData.Lz / max(obj.meshData.Nz - 1, 1);
            if dy == 0
                dy = 1;
                dz = 1;
            elseif dz == 0
                dz = 1;
            end
            fprintf(fileId, "# vtk DataFile Version 3.0\n");
            fprintf(fileId, "%s\n", obj.headerName);
            fprintf(fileId, "ASCII\n");
            fprintf(fileId, "DATASET STRUCTURED_POINTS\n");
            fprintf(fileId, "DIMENSIONS %d %d %d\n", obj.meshData.Nx, obj.meshData.Ny, obj.meshData.Nz);
            fprintf(fileId, "SPACING %f %f %f\n", dx, dy, dz);
            fprintf(fileId, "ORIGIN %d %d %d\n", 0, 0, 0);
        end
        
        function writeGmshHeader(obj, fileId)
            globalCoordinates = getGlobalCoordinates(obj.meshData);
            connectivity = getElements(obj.meshData);
            % Subtract 1 from the connectivity as VTK index from 0
            cells = [size(connectivity, 1)*ones(size(connectivity, 2), 1) connectivity'-1];
            
            fprintf(fileId, "# vtk DataFile Version 3.0\n");
            fprintf(fileId, "%s\n", obj.headerName);
            fprintf(fileId, "ASCII\n");
            fprintf(fileId, "DATASET UNSTRUCTURED_GRID\n");
            fprintf(fileId, "POINTS %d float\n", size(globalCoordinates, 2));
            fprintf(fileId, "%f %f %f\n", globalCoordinates(1, :), ...
                globalCoordinates(2, :), globalCoordinates(3, :));
            fprintf(fileId, "CELLS %d %d\n", size(cells, 1), numel(cells));
            fprintf(fileId, [repmat('%d ', size(cells)), "\n"], cells);
            fprintf(fileId, "CELL_TYPES %d\n", size(cells, 1));
            fprintf(fileId, "%d\n", obj.elementType.vtkCellType*ones(size(cells, 1), 1));
        end
        
        function writeDataset(obj, fileId, fieldData)
            switch(fieldData.type)
                case "scalars"
                    fprintf(fileId, "SCALARS %s float %d\n", ...
                        fieldData.name, fieldData.nbrComponents);
                    fprintf(fileId, "LOOKUP_TABLE default\n");
                case "vectors"
                    fprintf(fileId, "VECTORS %s float\n", fieldData.name);
            end
            fprintf(fileId, "%.8f\n", fieldData.data(:));
        end
    end
end

