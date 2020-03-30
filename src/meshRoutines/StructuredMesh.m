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
        Dimensions
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
            switch obj.Dimensions
                case 3
                    type = Elements.HEX_8;
                case 2
                    type = Elements.QUA_4;
                case 1
                    type = Elements.LIN_2;
            end
        end
        
        function dim = get.Dimensions(obj)
            if obj.Ny > 1 && obj.Nz > 1
                dim = 3;
            elseif obj.Ny > 1
                dim = 2;
            else
                dim = 1;
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
            [ix, iy, iz] = obj.toCartesian(elementNbr(:), true);
            
            hx = obj.Lx / (obj.Nx-1);
            dx = ceil(radius / hx);
            ixRange = max(1, ix - dx):min(obj.Nx-1, ix + dx);
            
            if obj.Ny == 1
                iyRange = 1;
            else    
                hy = obj.Ly / (obj.Ny-1);
                dy = ceil(radius / hy);
                iyRange = max(1, iy - dy):min(obj.Ny-1, iy + dy);
            end
            if obj.Nz == 1
                izRange = 1;
            else
                hz = obj.Lz / (obj.Nz-1);
                dz = ceil(radius / hz);
                izRange = max(1, iz - dz):min(obj.Nz-1, iz + dz);
            end

            neighbors = zeros(length(ixRange)*length(iyRange)*length(izRange), 1);
            c = 1;
            for iix = ixRange
                for iiy = iyRange
                    for iiz = izRange
                        neighbors(c) = obj.toLinear(iix, iiy, iiz, true);
                        c = c + 1;
                    end
                end
            end
        end

        function [n] = toLinear(obj, ix, iy, iz, element)
            if nargin == 5 && element
                Nx = max(obj.Nx - 1, 1);
                Ny = max(obj.Ny - 1, 1);
                Nz = max(obj.Nz - 1, 1);
            else
                Nx = obj.Nx;
                Ny = obj.Ny;
                Nz = obj.Nz;
            end
            if nargin < 3
                % Allow for the input to be 1 matrix with multiple columns
                if size(ix, 2) > 1
                    iy = ix(:, 2);
                    iz = ix(:, 3);
                    ix = ix(:, 1);
                else
                    iy = 1*ones(size(ix));
                    iz = 1*ones(size(ix));
                end
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

            
            n = (iz-1) * (Nx * Ny) + ...
                (iy-1) * (Nx) + ...
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
            
            if nargout == 1
                ix = [ix, iy, iz];
            end
        end
        
        function weights = elementWeights(obj, radius, weightFunction)
            middleElement = obj.toLinear(...
                max(ceil(obj.Nx-1)/2, 1), ...
                max(ceil(obj.Ny-1)/2, 1), ...
                max(ceil(obj.Nz-1)/2, 1), true);
            middleNeighbors = obj.elementNeighbors(middleElement, radius);
            maxNeighbors = numel(middleNeighbors);
            
            ne = obj.NumElements;
            I = zeros(ne*maxNeighbors, 1);
            J = zeros(ne*maxNeighbors, 1);
            V = zeros(ne*maxNeighbors, 1);
            % Calculate the weights
            c = 0;
            for e1 = 1:ne
                % Get the lower right node of the element to compute
                % distances
                elementNode = obj.toLinear(obj.toCartesian(e1, true));
                coord1 = obj.coordinates(elementNode);
                neighbors = obj.elementNeighbors(e1, radius);
                % Get the lower right nodes of the neighbor to compute
                % distances
                neighborNodes = obj.toLinear(obj.toCartesian(neighbors, true));
                coord2s = obj.coordinates(neighborNodes);
                dx = coord1(1) - coord2s(1, :);
                dy = coord1(2) - coord2s(2, :);
                dz = coord1(3) - coord2s(3, :);
                weights = weightFunction(dx, dy, dz);
                mask = weights > 0;
                indices = c + (1:nnz(mask));
                I(indices) = e1;
                J(indices) = neighbors(mask);
                V(indices) = weights(mask);
                c = c + nnz(mask);
            end
            weights = sparse(I(1:c), J(1:c), V(1:c), ne, ne);
        end
    end
end

