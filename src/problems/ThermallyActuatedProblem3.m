classdef ThermallyActuatedProblem3 < TopOptProblem
    %THERMALLYACTUATEDPROBLEM Summary of this class goes here
    %   Detailed explanation goes here
    
    methods
        function obj = ThermallyActuatedProblem3(flexFemModel, stiffFemModel, options, volumeFraction, u_max, weight, intermediateFunc)
            %THERMALLYACTUATEDPROBLEM Construct an instance of this class
            %   Detailed explanation goes here
            if nargin < 7
                intermediateFunc = [];
            end
            if nargin < 6
                weight = 0;
            end
            obj = obj@TopOptProblem(flexFemModel, options, intermediateFunc);
            stiffFemModel.assemble();
            obj.options.stiffFemModel = stiffFemModel;
            obj.options.u_max = u_max;
            obj.options.volumeFraction = volumeFraction;
            obj.options.intermediateWeight = weight;
        end
        
        function g = objective(obj, designPar)
            designPar = reshape(designPar, [], 1);
            %METHOD1 Summary of this method goes here
            %   Detailed explanation goes here
            filteredPar = obj.filterParameters(designPar);
            
            I = obj.fem.mechFEM.getDummy("josse");
            
            obj.fem.reassemble(filteredPar);
            obj.fem.solve();
            
            g = -sum(dot(I, obj.fem.mechFEM.displacements)) + ...
                obj.options.intermediateWeight*sum(filteredPar.*(1-filteredPar)) / numel(filteredPar);
        end
        
        function dgdphi = gradObjective(obj, designPar)
            designPar = reshape(designPar, [], 1);
            filteredPar = obj.filterParameters(designPar);

            I = obj.fem.mechFEM.getDummy("josse");
            
            adjointLoads_therm = zeros(size(obj.fem.temperatures));
            adjointLoads_mech = -I;
            
            dgdphi = obj.fem.gradChainTerm(adjointLoads_therm, adjointLoads_mech) + ...
                obj.options.intermediateWeight*(1-2*filteredPar) / numel(filteredPar);
            
            dgdphi(:) = obj.filterGradient(dgdphi, designPar);
        end
        
        function gs = constraint1(obj, designPar)
            designPar = reshape(designPar, [], 1);
            filteredPar = obj.filterParameters(designPar);
            % Stiff FEM model can be a mech FEM model
            st = obj.options.stiffFemModel;
            st.reassemble(filteredPar);
            st.solve();
            
            gs = st.displacements' * st.K_tot * st.displacements / obj.options.u_max - 1;
        end
        
        function dgsdphi = gradConstraint1(obj, designPar)
            designPar = reshape(designPar, [], 1);
            
            filteredPar = obj.filterParameters(designPar);
            
            dgsdphi = obj.gradCompliance(filteredPar) / obj.options.u_max;

            dgsdphi(:) = obj.filterGradient(dgsdphi, designPar);
        end
        
        function gs = constraint2(obj, designPar)
            designPar = reshape(designPar, [], 1);
            filteredPar = obj.filterParameters(designPar);
            % Stiff FEM model can be a mech FEM model
            gs = dot(filteredPar, obj.fem.volumes) / (obj.options.volumeFraction*sum(obj.fem.volumes)) - 1;
        end
        
        function dgsdphi = gradConstraint2(obj, designPar)
            designPar = reshape(designPar, [], 1);
            %filteredPar = obj.filterParameters(designPar);
            dgsdphi = obj.fem.volumes / (obj.options.volumeFraction * sum(obj.fem.volumes));

            dgsdphi(:) = obj.filterGradient(dgsdphi, designPar);
        end
    end
    
    methods(Access = protected)
        function dgdphi = gradCompliance(obj, designPar)
            designPar = reshape(designPar, [], 1);

            st = obj.options.stiffFemModel;
            adjointLoads = 2 * st.K_tot * st.displacements;
            
            dgdphi = st.gradChainTerm(adjointLoads);
            
            for e = 1:length(designPar)
                d = designPar(e);
                dEdphi = st.stiffnessDer(d);
                k0 = st.getElementBaseMatrix(e, 'D');
                u_e = st.displacements(st.Edof(:, e), :);
                kT = dEdphi * k0*u_e;
                
                dgdphi(e) = dgdphi(e) + sum(dot(u_e, kT));
            end
        end
    end
end

