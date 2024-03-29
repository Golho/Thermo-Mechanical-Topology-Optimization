classdef ThermallyActuatedProblemTransient < TopOptProblem
    %THERMALLYACTUATEDPROBLEM Summary of this class goes here
    %   Detailed explanation goes here
    
    methods
        function obj = ThermallyActuatedProblemTransient(flexFemModel, stiffFemModel, ...
                dummyName, options, massLimit, u_max, intermediateFunc)
            %THERMALLYACTUATEDPROBLEM Construct an instance of this class
            %   Detailed explanation goes here
            if nargin < 7
                intermediateFunc = [];
            end
            obj = obj@TopOptProblem(flexFemModel, options, intermediateFunc);
            stiffFemModel.assemble();
            obj.options.stiffFemModel = stiffFemModel;
            obj.options.dummyName = dummyName;
            obj.options.massLimit = massLimit;
            obj.options.u_max = u_max;
            [obj.options.densityFunc, obj.options.densityDerFunc] = ...
                densitySIMP(obj.options.materials, 1);
        end
        
        function error = objective(obj, designPar)
            %METHOD1 Summary of this method goes here
            %   Detailed explanation goes here
            obj.fem.reassemble(designPar);
            obj.fem.solve();
            
            I = obj.fem.mechFEM.getDummy(obj.options.dummyName);
            i = logical(I);
            diff = obj.fem.mechFEM.displacements(i) - I(i);
            error = sum(diff.^2, "all") / obj.options.u_max^2;
        end
        
        function dgdphi = gradObjective(obj, designPar)
            %filteredPar = obj.filterParameters(designPar);

            I = obj.fem.mechFEM.getDummy(obj.options.dummyName);
            
            adjointLoads_therm = zeros(size(obj.fem.temperatures));
            diff = zeros(size(I));
            i = logical(I);
            diff(i) = obj.fem.mechFEM.displacements(i) - I(i);
            adjointLoads_mech = 2*diff / obj.options.u_max^2;
            
            dgdphi = obj.fem.gradChainTerm(adjointLoads_therm, adjointLoads_mech);
        end
        
        function gs = constraint0(obj, designPar)
            % Stiff FEM model can be a mech FEM model
            st = obj.options.stiffFemModel;
            st.reassemble(designPar);
            st.solve();
            
            gs = st.displacements' * st.K_tot * st.displacements / obj.options.u_max - 1;
        end
        
        function dgsdphi = gradConstraint0(obj, designPar)
            dgsdphi = obj.gradCompliance(designPar) / obj.options.u_max;
        end
        
        function gs = constraint1(obj, designPar)
            densities = obj.options.densityFunc(designPar);
            gs = dot(densities, obj.fem.volumes) / obj.options.massLimit - 1;
        end
        
        function dgsdphi = gradConstraint1(obj, designPar)
            dgsdphi = obj.options.densityDerFunc(designPar) .* obj.fem.volumes' / obj.options.massLimit;
        end
        
        function plotResults(obj)
            % Get the nodes which to plot the path of
            I = obj.fem.mechFEM.getDummy(obj.options.dummyName);
            controlDofs = any(I, 2);
            controlNodes = find(any(reshape(controlDofs, obj.fem.mechFEM.spatialDimensions, []), 1));
            
            for controlNode = controlNodes
                figure
                hold on
                pathPlot(obj.fem.mechFEM.getDofs(controlNode), obj.fem.mechFEM.displacements)
                pathPlot(obj.fem.mechFEM.getDofs(controlNode), I, '+');
                hold off
            end
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

