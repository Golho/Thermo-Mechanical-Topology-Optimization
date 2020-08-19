classdef ThermallyActuatedProblem3_Robust < RobustProblem
    %THERMALLYACTUATEDPROBLEM Summary of this class goes here
    %   Detailed explanation goes here
    
    properties(Access = protected)
        maxC2Index
    end
    
    methods
        function obj = ThermallyActuatedProblem3_Robust(flexFemModel, stiffFemModel, options, massLimit, u_max, filters, intermediateFunc)
            %THERMALLYACTUATEDPROBLEM Construct an instance of this class
            %   Detailed explanation goes here
            if nargin < 7
                intermediateFunc = [];
            end
            obj = obj@RobustProblem(flexFemModel, options, massLimit, filters, intermediateFunc);
            stiffFemModel.assemble();
            obj.options.stiffFemModel = stiffFemModel;
            obj.options.u_max = u_max;
            [obj.options.densityFunc, obj.options.densityDerFunc] = ...
                densitySIMP(obj.options.materials, ones(size(obj.fem.designPar, 1), 1));
        end
        
        function g = subObjective(obj, designPar)
            %METHOD1 Summary of this method goes here
            %   Detailed explanation goes here
            
            I = obj.fem.mechFEM.getDummy("josse");
            
            obj.fem.reassemble(designPar);
            obj.fem.solve();
            
            g = -sum(dot(I, obj.fem.mechFEM.displacements));
        end
        
        function dgdphi = subGradObjective(obj, designPar)

            I = obj.fem.mechFEM.getDummy("josse");
            
            adjointLoads_therm = zeros(size(obj.fem.temperatures));
            adjointLoads_mech = -I;
            
            dgdphi = obj.fem.gradChainTerm(adjointLoads_therm, adjointLoads_mech);
        end
        
        function gs = volumeConstraint(obj, designPar)
            densities = obj.options.densityFunc(designPar);
            gs = dot(densities, obj.fem.volumes) / obj.dilatedMassLimit - 1;
        end
        
        function dgsdphi = volumeGradConstraint(obj, designPar)
            dgsdphi = obj.options.densityDerFunc(designPar) .* obj.fem.volumes' / obj.dilatedMassLimit;
        end
        
        function gs = constraint2(obj, designPar)
            thresholdDesigns = cell(1, 3);
            gss = zeros(3, 1);
            for i = 1:3
                thresholdDesigns{i} = obj.options.filters(i).filter(designPar);
                gss(i) = obj.subConstraint2(thresholdDesigns{i});
            end
            
            [gs, obj.maxC2Index] = max(gss);
        end
        
        function gs = subConstraint2(obj, designPar)
            % Stiff FEM model can be a mech FEM model
            st = obj.options.stiffFemModel;
            st.reassemble(designPar);
            st.solve();
            
            gs = st.displacements' * st.K_tot * st.displacements / obj.options.u_max - 1;
        end
        
        function dgsdphi = gradConstraint2(obj, designPar)
            projectedPar = ...
                obj.options.filters(obj.maxC2Index).filter(designPar);
            dgsdphi = obj.gradCompliance(projectedPar) / obj.options.u_max ...
                .* obj.options.filters(obj.maxC2Index).gradFilter(designPar);
        end
    end
    
    methods(Access = protected)
        function dgdphi = gradCompliance(obj, designPar)
            st = obj.options.stiffFemModel;
            adjointLoads = 2 * st.K_tot * st.displacements;
            
            dgdphi = st.gradChainTerm(adjointLoads);
            
            for e = 1:size(designPar, 2)
                d = designPar(:, e);
                dEdphi = st.stiffnessDer(d);
                k0 = st.getElementBaseMatrix(e, 'D');
                u_e = st.displacements(st.Edof(:, e), :);
                kT = k0*u_e;
                
                dgdphi(:, e) = dgdphi(:, e) + dEdphi * sum(dot(u_e, kT));
            end
        end
    end
end

