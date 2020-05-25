classdef OptMechFEMStructured < MechFEMStructured & OptMechFEMBase
    %OPTMECHFEMSTRUCTURED Summary of this class goes here
    %   Detailed explanation goes here
    
    properties
        volumes
    end
    
    methods
        function obj = OptMechFEMStructured(varargin)
            %OPTMECHFEMSTRUCTURED Construct an instance of this class
            %   Detailed explanation goes here

            obj = obj@MechFEMStructured(varargin{:});
            obj.designPar = ones(obj.mesh.NumElements, 1);
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
            
            k0_uu = obj.getElementBaseMatrix(1, "D");
            obj.K = obj.optIntegrate(k0_uu, "E", 2);
            
            k0_uT = obj.getElementBaseMatrix(1, "D-alpha");
            obj.K_thermal = obj.optIntegrate(k0_uT, "alpha", 3);
            
            obj.assemble();
        end
        
        function weights = computeWeights(obj, radius, weightFunction)
            startTime = tic;
            
            weights = obj.mesh.elementWeights(radius, weightFunction);
            
            % Normalize so the sum of the weights equal 1
            weights = weights ./ sum(weights, 2);
            fprintf("Computed weights:\t\t\t%f secs\n", toc(startTime));
        end
        
        function chainGrad = gradChainTerm(obj, adjointLoads)
            % Get chain rule term for the sensitivity for a mechanical
            % system where the temperature changes are not dependent on the
            % design parameter (otherwise a coupled system must be used)
            startTime = tic;
            
            adjoints = obj.solveAdjoint(adjointLoads);
            
            chainGrad = zeros(size(obj.designPar));
            
            k0_uu = obj.getElementBaseMatrix(1, 'D');
            k0_uT = obj.getElementBaseMatrix(1, 'D-alpha');
            
            edof = obj.Edof(:, 1);
            enod = obj.Enod(:, 1);
            u_e = obj.displacements(edof, :);
            T_e = obj.temperatureChanges(enod, :);
            adjoint_e = adjoints(edof, :);
            dRdx = zeros(size(adjoint_e));
            
            for e = 1:length(chainGrad)
                edof(:) = obj.Edof(:, e);
                enod(:) = obj.Enod(:, e);
                u_e(:) = obj.displacements(edof, :);
                T_e(:) = obj.temperatureChanges(enod, :);
                adjoint_e(:) = adjoints(edof, :);
            
                dEdphi = obj.stiffnessDer(obj.designPar(e));
                dalphadphi = obj.thermalExpDer(obj.designPar(e));

                dRdx(:) = dEdphi * k0_uu * u_e - dalphadphi * k0_uT * T_e;
                chainGrad(e) = -sum(dot(adjoint_e, dRdx));
            end
            fprintf("Computed gradient term:\t\t%f secs\n", toc(startTime));
        end
    end
    
    
    methods(Access = protected)
        function [matrix] = optIntegrate(obj, elementMatrix, property, dim)
            nbrElements = obj.mesh.NumElements;
            
            switch dim
                case 2
                    EdofEdof = repmat(obj.Edof, size(obj.Edof, 1), 1);
                    I = EdofEdof(:);
                    EEddooff = repelem(obj.Edof(:), size(obj.Edof, 1), 1);
                    J = EEddooff;
                    rows = obj.nbrDofs;
                    cols = obj.nbrDofs;
                case 3
                    EdofEdof = repmat(obj.Edof, size(obj.Enod, 1), 1);
                    I = EdofEdof(:);
                    EEnnoodd = repelem(obj.Enod(:), size(obj.Edof, 1), 1);
                    J = EEnnoodd;
                    rows = obj.nbrDofs;
                    cols = obj.nbrNodes;
                otherwise
                    error("Dim must be either 2 or 3");
            end
            
            V = zeros(size(I));
            c = 0; % offset
            for e = 1:nbrElements
                d = obj.designPar(e);
                factor = obj.elementProp(property, d);
                V(c+(1:numel(elementMatrix))) = factor * ...
                    elementMatrix(:);
                c = c + numel(elementMatrix);
            end
            matrix = sparse(I, J, V, rows, cols);
        end
    end
end

