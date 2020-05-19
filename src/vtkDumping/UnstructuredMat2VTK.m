function UnstructuredMat2VTK(filename, coords, connectivity, elementType, ...
    dataMatrix, format, label, header, cellData)
% Writes a matrix as a *.VTK file as input for Paraview.

assert(size(coords, 1) == 3, "The coordinate matrix must be a [3 x N] matrix");

if nargin < 9
    cellData = false;
end

[Nx, Ny, Nz] = size(dataMatrix);
% Subtract 1 from the connectivity as VTK index from 0
cells = [size(connectivity, 1)*ones(size(connectivity, 2), 1) connectivity'-1];

% Open the file.
fid = fopen(filename, "w", "n", "US-ASCII");
if fid == -1
    error("Cannot open file for writing.");
end

switch format
    case "ascii"
        fprintf(fid, "# vtk DataFile Version 3.0\n");
        fprintf(fid, "%s\n", header);
        fprintf(fid, "ASCII\n");
        fprintf(fid, "DATASET UNSTRUCTURED_GRID\n");
        fprintf(fid, "POINTS %d float\n", size(coords, 2));
        fprintf(fid, "%f %f %f\n", coords);
        fprintf(fid, "CELLS %d %d\n", size(cells, 1), numel(cells));
        fprintf(fid, [repmat('%d ', 1, size(cells, 2)), '\n'], cells');
        fprintf(fid, "CELL_TYPES %d\n", size(cells, 1));
        fprintf(fid, "%d\n", elementType.vtkCellType*ones(size(cells, 1), 1));
        if cellData
            fprintf(fid, "CELL_DATA %d\n", Nx*Ny*Nz);
        else
            fprintf(fid, "POINT_DATA %d\n", Nx*Ny*Nz);
        end
        fprintf(fid, "SCALARS %s float 1\n", label);
        fprintf(fid, "LOOKUP_TABLE default\n");
        fprintf(fid, "%f\n", dataMatrix(:));
    case "binary"
        fprintf(fid, "# vtk DataFile Version 3.0\n");
        fprintf(fid, "%s\n", header);
        fprintf(fid, "ASCII\n");
        fprintf(fid, "DATASET UNSTRUCTURED_GRID\n");
        fprintf(fid, "POINTS %d float\n", size(coords, 2));
        fprintf(fid, "%f %f %f\n", coords(1, :), coords(2, :), coords(3, :));
        fprintf(fid, "CELLS %d %d\n", size(cells, 1), numel(cells));
        fprintf(fid, [repmat('%d ', size(cells)), "\n"], cells);
        fprintf(fid, "CELL_TYPES %d\n", size(cells, 1));
        fprintf(fid, "%d\n", elementType.vtkCellType*ones(size(cells, 1), 1));
        if cellData
            fprintf(fid, "DIMENSIONS %d %d %d\n", Nx+1, Ny+1, Nz+1);
        else
            fprintf(fid, "DIMENSIONS %d %d %d\n", Nx, Ny, Nz);
        end
        fprintf(fid, "SPACING %d %d %d\n", 1, 1, 1);
        fprintf(fid, "ORIGIN %d %d %d\n", 0, 0, 0);
        if cellData
            fprintf(fid, "CELL_DATA %d\n", Nx*Ny*Nz);
        else
            fprintf(fid, "POINT_DATA %d\n", Nx*Ny*Nz);
        end
        fprintf(fid, "SCALARS %s float 1\n", label);
        fprintf(fid, "LOOKUP_TABLE default\n");
        fwrite(fid, dataMatrix(:), "float", "ieee-be");
    otherwise
        error("The format argument must either be 'ascii' or 'binary'");
end

% Close the file.
fclose(fid);
end
