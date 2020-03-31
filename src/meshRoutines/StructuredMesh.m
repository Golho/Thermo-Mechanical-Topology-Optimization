classdef StructuredMesh
    %UNTITLED Summary of this class goes here
    %   Detailed explanation goes here
    
    properties
        Nx
        Ny
        Nz
        Lx
        Ly
        Lz
    end
    
    properties (Dependent)
        NumNodes
        NumElements
        ElementType
    end
    
    methods
        function obj = StructuredMesh(xData, yData, zData)
            %UNTITLED Construct an instance of this class
            %   Detailed explanation goes here
            obj.Nx = xData(1);
            obj.Lx = xData(2);
            
            if exist('yData', 'var')
                obj.Ny = yData(1);
                obj.Ly = yData(2);
            else
                obj.Ny = 1;
                obj.Ly = 0;
            end
            
            if exist('zData', 'var')
                obj.Nz = zData(1);
                obj.Lz = zData(2);
            else
                obj.Nz = 1;
                obj.Lz = 0;
            end
            
            fprintf("Created a structured mesh (%d nodes, %d elements)\n", ...
                obj.NumNodes, obj.NumElements);
        end
        
        function numE = get.NumElements(obj)
            %METHOD1 Summary of this method goes here
            %   Detailed explanation goes here
            numE = (obj.Nx-1)*max(1, obj.Ny-1)*max(1, obj.Nz-1);
        end
        
        function numN = get.NumNodes(obj)
            numN = obj.Nx*obj.Ny*obj.Nz;
        end
        
        function type = get.ElementType(obj)
            if obj.Ny > 1 && obj.Nz > 1
                type = Elements.HEX_8;
            elseif obj.Ny > 1
                type = Elements.QUA_4;
            else
                type = Elements.LIN_2;
            end
        end
        
        function coord = coordinates(obj, n)
            if nargin == 1
                n = 1:obj.NumNodes;
            end
            [ix, iy, iz] = obj.toCartesian(n(:));
            coord = zeros(3, length(n));
            coord(1, :) = (ix-1) * obj.Lx / (obj.Nx-1);
            coord(2, :) = (iy-1) * obj.Ly / (max(obj.Ny-1, 1));
            coord(3, :) = (iz-1) * obj.Lz / (max(obj.Nz-1, 1));
        end
        
        function [Ex, Ey, Ez] = elementCoordinates(obj)
            connectivity = obj.elementNodes(1:obj.NumElements);
            coordinates = obj.coordinates();
            x = coordinates(1, :);
            y = coordinates(2, :);
            z = coordinates(3, :);
            Ex = x(connectivity);
            Ey = y(connectivity);
            Ez = z(connectivity);
        end
        
        function connectivity = elementNodes(obj, elementNbr)
            assert(all(elementNbr >= 1 & elementNbr <= obj.NumElements), ...
                "Element number must be within range");
            
            % Imagine every node "owning" an element, except the ones at
            % the end of a line, i.e. when ix == obj.Nx, iy == obj.Ny etc.
            [ix, iy, iz] = obj.toCartesian(elementNbr(:), true);
            if obj.Nz == 1 && obj.Ny == 1
                ixs = [ix, ix+1];
                iys = [iy, iy];
                izs = [iz, iz];
                connectivity = obj.toLinear(ixs, iys, izs)';
            elseif obj.Nz == 1
                ixs = [ix, ix+1, ix+1, ix];
                iys = [iy, iy, iy+1, iy+1];
                izs = [iz, iz, iz, iz];
                connectivity = obj.toLinear(ixs, iys, izs)';
            else
                ixs = [ix, ix+1, ix+1, ix, ix, ix+1, ix+1, ix];
                iys = [iy, iy, iy+1, iy+1, iy, iy, iy+1, iy+1];
                izs = [iz, iz, iz, iz, iz+1, iz+1, iz+1, iz+1];
                connectivity = obj.toLinear(ixs, iys, izs)';
            end
        end
        
        function [neighbors] = elementNeighbors(obj, elementNbr, radius)
            hx = obj.Lx / (obj.Nx-1);
            dx = ceil(radius / hx);
            
            if obj.Ny == 1
                dy = 0;
            else    
                hy = obj.Ly / (obj.Ny-1);
                dy = ceil(radius / hy);
            end
            if obj.Nz == 1
                dz = 0;
            else
                hz = obj.Lz / (obj.Nz-1);
                dz = ceil(radius / hz);
            end
            
            

            [ix, iy, iz] = obj.toCartesian(elementNbr(:), true);
            ixRange = max(1, ix - dx):min(obj.Nx, ix + dx);
            iyRange = max(1, iy - dy):min(obj.Ny, iy + dy);
            izRange = max(1, iz - dz):min(obj.Nz, iz + dz);
            neighbors = zeros(length(ixRange)*length(iyRange)*length(izRange), 1);
            c = 1;
            for iix = ixRange
                for iiy = iyRange
                    for iiz = izRange
                        neighbors(c) = obj.toLinear(iix, iiy, iiz);
                        c = c + 1;
                    end
                end
            end
        end

        function [n] = toLinear(obj, ix, iy, iz)
            if nargin < 3
                iy = 1*ones(size(ix));
                iz = 1*ones(size(ix));
            elseif nargin < 4
                iz = 1*ones(size(ix));
            end
            % Comment out the assertions, as they are time-consuming
%             assert(all(mod([ix, iy, iz], 1) == 0, 'all'), ...
%                 "Cartesian indices must be integers");
%             xInRange = all(ix > 0 & ix <= obj.Nx, 'all');
%             yInRange = all(iy > 0 & iy <= obj.Ny, 'all');
%             zInRange = all(iz > 0 & iz <= obj.Nz, 'all');
%             assert(xInRange && yInRange && zInRange, ...
%                 "Cartesian indices must be within range");
            
            n = (iz-1) * (obj.Nx * obj.Ny) + ...
                (iy-1) * (obj.Nx) + ...
                ix;
        end
        
        function [ix, iy, iz] = toCartesian(obj, n, element)
            if nargin == 3 && element
                Nx = max(obj.Nx - 1, 1);
                Ny = max(obj.Ny - 1, 1);
                Nz = max(obj.Nz - 1, 1);
            else
                Nx = obj.Nx;
                Ny = obj.Ny;
                Nz = obj.Nz;
            end
            % Comment out the assertions, as they are time-consuming
%             assert(all(n >= 1 & n <= obj.NumNodes, 'all'), "Index must be within range");
            iz = ceil(n ./ (Nx * Ny));
            iy = mod(ceil(n ./ Nx)-1, Ny) + 1;
            ix = mod(n-1, Nx)+1;
        end
    end
end

