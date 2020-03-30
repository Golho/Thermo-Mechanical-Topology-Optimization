function Mat2VTK(filename, matrix, format, label, header, cellData)
% Writes a matrix as a *.VTK file as input for Paraview.

if nargin < 6
    cellData = false;
end

[Nx, Ny, Nz] = size(matrix);

% Open the file.
fid = fopen(filename, "w");
if fid == -1
    error("Cannot open file for writing.");
end

switch format
    case "ascii"
        fprintf(fid, "# vtk DataFile Version 3.0\n");
        fprintf(fid, "%s\n", header);
        fprintf(fid, "ASCII\n");
        fprintf(fid, "DATASET STRUCTURED_POINTS\n");
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
        fprintf(fid, "%f\n", matrix(:));
    case "binary"
        fprintf(fid, "# vtk DataFile Version 3.0\n");
        fprintf(fid, "%s\n", header);
        fprintf(fid, "BINARY\n");
        fprintf(fid, "DATASET STRUCTURED_POINTS\n");
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
        fwrite(fid, matrix(:), "float", "ieee-be");
    otherwise
        error("The format argument must either be 'ascii' or 'binary'");
end

% Close the file.
fclose(fid);
end
