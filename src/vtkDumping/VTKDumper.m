classdef VTKDumper < handle
    %VTKDumper Class to dump data to .vtk-files
    %   This class allows for adding several data sets on a Gmsh-mesh or 
    %   StructuredMesh and then dump everything into a set of .VTK-files 
    %   (one for each time step)
    
    properties
        headerName              % Name of the header in the .VTK
        meshData                % Mesh (Gmsh or StructuredMesh)
        structured              % Boolean if mesh is structured
        elementType             % Element type (for unstructured meshes)
        pointData = struct(...  % Array of point data on the mesh
            "name", [], ...
            "type", [], ...
            "data", [], ...
            "nbrComponents", [], ...
            "timeSteps", [] ...
            );
        nbrPointData = 0;
        cellData = struct(...   % Array of cell data on the mesh
            "name", [], ...
            "type", [], ...
            "data", [], ...
            "nbrComponents", [], ...
            "timeSteps", [] ...
            );
        nbrCellData = 0;
        timeSteps = 1;          % Number of time steps
        
        extrudeZ = false;       % Boolean for if a 2D mesh should be extruded into 3D
    end
    
    methods
        function obj = VTKDumper(headerName, meshData, elementType)
            if isa(meshData, "Gmsh")
                assert(nargin == 3, "The element type must be specified for a GMSH mesh");
                obj.structured = false;
                obj.elementType = elementType;
            elseif isa(meshData, "StructuredMesh")
                assert(nargin ~= 3, "The element type is superfluous when mesh is structured");
                obj.structured = true;
            else
                error("The identified mesh type is not supported");
            end
            obj.headerName = headerName;
            obj.meshData = meshData;
        end
        
        function addData(obj, fieldData, dataType)
            %ADDDATA Add point or cell data on the mesh
            %   addData(obj, fieldData, dataType) Add fieldData with
            %   dataType ("point" or "cell") to dump later. The structure
            %   of the fieldData can been in the properties of the
            %   VTKDumper.
            
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
            % Update the number of time steps
            obj.timeSteps = max(max(fieldData.timeSteps), obj.timeSteps);
        end
        
        function dump(obj, filePrefix)
            %DUMP   Dump the data added to the VTKDumper into a set of
            %.VTK-files
            %   dump(obj, filePrefix) Dump the data to .VTK-files with the
            %   prefix filePrefix and a suffix showing the time step. E.g.
            %   "stress_4.vtk"
            for timeStep = 1:obj.timeSteps
                filename = filePrefix + "_" + timeStep + ".vtk";
                % Open the file.
                fid = fopen(filename, "w");
                if fid == -1
                    error("Cannot open file for writing.");
                end
                
                % Write the header
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
                
                % Write all the cell data
                if numel(obj.cellData) > 0
                    fprintf(fid, "CELL_DATA %d\n", obj.meshData.NumElements);
                    for fieldData = obj.cellData
                        if ismember(timeStep, fieldData.timeSteps)
                            obj.writeDataset(fid, fieldData);
                        end
                    end
                end
                
                % Write all the point data
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
            % Write the header for a .vtk for a structured mesh
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
            % Write the header for a .vtk for a unstructured mesh
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
            % Write the data of fieldData into a .vtk-file
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

