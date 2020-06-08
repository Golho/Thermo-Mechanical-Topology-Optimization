classdef OptHeatFEM < HeatFEM & OptHeatFEMBase
    %HEATFEM Summary of this class goes here
    %   Detailed explanation goes here
    
    properties
        volumes
    end
    
    methods
        function obj = OptHeatFEM(nbrMaterials, varargin)
            %HEATFEM Construct an instance of this class
            %   Detailed explanation goes here
            obj = obj@HeatFEM(varargin{:});
            obj.designPar = ones(nbrMaterials-1, size(obj.Enod, 2));
        end
        
        function volumes = get.volumes(obj)
            volumes = obj.elementVolume(obj.Ex, obj.Ey, obj.Ez, obj.ElementType)';
        end
        
        function reassemble(obj, designPar)
            obj.designPar = designPar;
            [obj.K, obj.M] = obj.applyMain();

            if obj.transient
                dt = obj.tFinal / (obj.timeSteps-1);
                obj.A = obj.M + dt*obj.theta*obj.K;
                obj.B = obj.M + dt*(obj.theta - 1)*obj.K;
            end
        end
        
        function weights = computeWeights(obj, radius, weightFunction)
            startTime = tic;
            % Calculate the center point of each element
            nbrMainElems = size(obj.Enod, 2);
            I = zeros(nbrMainElems*10, 1);
            J = zeros(nbrMainElems*10, 1);
            V = zeros(nbrMainElems*10, 1);
            elementCoord = zeros(3, nbrMainElems);
            for e = 1:nbrMainElems
                enod = obj.Enod(:, e);
                elementCoord(:, e) = mean(obj.nodeCoordinates(:, enod), 2);
            end
            % Calculate the weights
            c = 0;
            for e1 = 1:nbrMainElems
                coord1 = elementCoord(:, e1);
                I(c+1) = e1;
                J(c+1) = e1;
                V(c+1) = weightFunction(0, 0, 0);
                c = c + 1;
                
                % Assume all "older" elements to be neighbors and then add
                % the symmetric counter-parts of all matches
                neighbors = 1:e1-1;
                coord2s = elementCoord(:, neighbors);
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
                indices = c + (1:nnz(mask));
                % Add the symmetric counter-part
                I(indices) = neighbors(mask);
                J(indices) = e1;
                V(indices) = weights(mask);
                c = c + nnz(mask);
            end
            weights = sparse(I(1:c), J(1:c), V(1:c));
            
            % Normalize so the sum of the weights equal 1
            weights = weights ./ sum(weights, 2);
            fprintf("Computed weights:\t\t\t%f secs\n", toc(startTime));
        end
        
        function chainGrad = gradChainTerm(obj, adjointLoads, dT0dx)
            startTime = tic;
            if nargin < 3
                dT0dx = sparse(size(obj.temperatures, 1), ...
                    numel(obj.designPar));
            end
            
            deltaT = obj.tFinal / (obj.timeSteps - 1);
            
            adjoints = obj.solveAdjoint(adjointLoads);
            
            chainGrad = (-adjoints(:, 1)' * obj.B * dT0dx)';
            chainGrad = reshape(chainGrad, size(obj.designPar));
            
            enod = obj.Enod(:, 1);
            T_e = obj.temperatures(enod, :);
            adjoint_e = adjoints(enod, :);
            adjoint_dRdx = zeros(size(obj.designPar, 1), 1);
            
            for e = 1:length(chainGrad)                
                enod(:) = obj.Enod(:, e);
                T_e(:) = obj.temperatures(enod, :);
                adjoint_e(:) = adjoints(enod, :);
            
                d = obj.designPar(:, e);
                dkappadphi = obj.conductivityDer(d);
                dcpdphi = obj.heatCapacityDer(d);
                
                k0 = obj.getElementBaseMatrix(1, 'D');
                c0 = obj.getElementBaseMatrix(1, 'cp');
                
                adjoint_dRdx(:) = dkappadphi * sum(dot(adjoint_e, ...
                        deltaT*obj.theta*k0*T_e(:, 2:end) - deltaT*(obj.theta-1)*k0*T_e(:, 1:end-1)...
                    )) + ...
                    dcpdphi * sum(dot(adjoint_e, c0 * (T_e(:, 2:end) - T_e(:, 1:end-1))));
                chainGrad(:, e) = chainGrad(:, e) - adjoint_dRdx;
            end
            fprintf("Computed gradient term:\t\t%f secs\n", toc(startTime));
        end
    end
    methods(Access = protected)
        function [matrix] = optIntegrate(obj, func, property)
            nbrElements = size(obj.Enod, 2);

            EnodEnod = repmat(obj.Enod, size(obj.Enod, 1), 1);
            I = EnodEnod(:);
            EEnnoodd = repelem(obj.Enod(:), size(obj.Enod, 1), 1);
            J = EEnnoodd;
            
            V = zeros(size(I));
            c = 0; % offset
            for e = 1:nbrElements
                d = obj.designPar(:, e);
                edof = obj.Enod(:, e);
                elementCoord = obj.nodeCoordinates(:, edof)';
                elementMatrix = func(elementCoord);
                V(c+(1:numel(elementMatrix))) = obj.elementProp(property, d) * ...
                    elementMatrix(:);
                c = c + numel(elementMatrix);
            end
            
            matrix = sparse(I, J, V, obj.nbrDofs, obj.nbrDofs);
        end
        
        function [K, M] = applyMain(obj)
            startTime = tic;
            stiffFunc = @(elementCoord) obj.elementStiffness(elementCoord, ...
                obj.ElementType, eye(3));
            K = obj.optIntegrate(stiffFunc, 'kappa');

            massFunc = @(elementCoord) obj.elementMass(elementCoord, ...
                obj.ElementType, 1);
            M = obj.optIntegrate(massFunc, 'cp');
            fprintf("Integrated main BC:\t\t\t%f secs\n", toc(startTime));
        end
    end
end
