classdef OptHeatFEMStructured < HeatFEMStructured & OptHeatFEMBase
    %HEATFEM Summary of this class goes here
    %   Detailed explanation goes here
    
    properties
        volumes
    end
    
    methods
        function obj = OptHeatFEMStructured(nbrMaterials, varargin)
            %HEATFEM Construct an instance of this class
            %   Detailed explanation goes here
            
            obj = obj@HeatFEMStructured(varargin{:});
            obj.designPar = ones(nbrMaterials-1, obj.mesh.NumElements);
        end
        
        function volumes = get.volumes(obj)
            ex = obj.Ex(:, 1);
            ey = obj.Ey(:, 1);
            ez = obj.Ez(:, 1);
            % Take advantage of the elements being identical
            volume = obj.elementVolume(ex, ey, ez, obj.ElementType);
            volumes = volume*ones(obj.mesh.NumElements, 1);
        end
        
        function reassemble(obj, designPar)
            obj.designPar = designPar;

            enod = obj.Enod(:, 1);
            elementCoord = obj.nodeCoordinates(:, enod)';
            [obj.K, obj.M] = obj.applyMain(elementCoord);

            if obj.transient
                dt = obj.tFinal / (obj.timeSteps-1);
                obj.A = obj.M + dt*obj.theta*obj.K;
                obj.B = obj.M + dt*(obj.theta - 1)*obj.K;
            end
        end
        
        function weights = computeWeights(obj, radius, weightFunction)
            startTime = tic;
            
            weights = obj.mesh.elementWeights(radius, weightFunction);
            
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
            
            k0 = obj.getElementBaseMatrix(1, 'D');
            c0 = obj.getElementBaseMatrix(1, 'cp');
            
            enod = obj.Enod(:, 1);
            T_e = obj.temperatures(enod, :);
            adjoint_e = adjoints(enod, :);
            adjoint_dRdx = zeros(size(obj.designPar, 1), 1);
            
            for e = 1:size(chainGrad, 2)
                d = obj.designPar(:, e);
                enod(:) = obj.Enod(:, e);
                T_e(:) = obj.temperatures(enod, :);
                adjoint_e(:) = adjoints(enod, :);
            
                dkappadphi = obj.conductivityDer(d);
                dcpdphi = obj.heatCapacityDer(d);
                
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
        function [matrix] = optIntegrate(obj, elementMatrix, property)
            nbrElements = obj.mesh.NumElements;

            EnodEnod = repmat(obj.Enod, size(obj.Enod, 1), 1);
            I = EnodEnod(:);
            EEnnoodd = repelem(obj.Enod(:), size(obj.Enod, 1), 1);
            J = EEnnoodd;
            
            V = zeros(size(I));
            c = 0; % offset
            for e = 1:nbrElements
                d = obj.designPar(:, e);
                V(c+(1:numel(elementMatrix))) = obj.elementProp(property, d) * ...
                    elementMatrix(:);
                c = c + numel(elementMatrix);
            end
            
            matrix = sparse(I, J, V, obj.nbrDofs, obj.nbrDofs);
        end
        
        function [K, M] = applyMain(obj, elementCoord)
            startTime = tic;
            k0 = obj.elementStiffness(elementCoord, ...
                obj.mesh.ElementType, eye(3));
            K = obj.optIntegrate(k0, 'kappa');

            c0 = obj.elementMass(elementCoord, obj.mesh.ElementType, 1);
            M = obj.optIntegrate(c0, 'cp');
            fprintf("Integrated main BC:\t\t\t%f secs\n", toc(startTime));
        end
    end
end
