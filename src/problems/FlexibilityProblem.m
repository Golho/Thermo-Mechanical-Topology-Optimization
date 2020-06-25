classdef FlexibilityProblem < TopOptProblem
    %FLEXIBILITYPROBLEM
    %   An dummy must be prescribed on the FEM model with the name "output"
    
    methods
        function obj = FlexibilityProblem(femModel, options, massLimit, u_max, intermediateFunc)
            %UNTITLED2 Construct an instance of this class
            %   Detailed explanation goes here
            if nargin < 5
                intermediateFunc = [];
            end
            obj = obj@TopOptProblem(femModel, options, intermediateFunc);
            obj.options.u_max = u_max;
            obj.options.massLimit = massLimit;
            [obj.options.densityFunc, obj.options.densityDerFunc] = densitySIMP(obj.options.materials, 1);
        end
        
        function g = objective(obj, designPar)
            %METHOD1 Summary of this method goes here
            %   Detailed explanation goes here
            I = obj.fem.getDummy("josse");
            
            obj.fem.reassemble(designPar);
            obj.fem.solve();

            g = - I' * obj.fem.displacements / obj.options.u_max;
        end
        
        function dgdphi = gradObjective(obj, designPar)
            
            %filteredPar = obj.filterParameters(designPar);

            I = obj.fem.getDummy("josse");
            
            adjointLoads = -I / obj.options.u_max;
            
            dgdphi = obj.fem.gradChainTerm(adjointLoads);
        end        
        
        function gs = constraint1(obj, designPar)
            densities = obj.options.densityFunc(designPar);
            gs = dot(densities, obj.fem.volumes) / obj.options.massLimit - 1;
        end
        
        function dgsdphi = gradConstraint1(obj, designPar)
            dgsdphi = obj.options.densityDerFunc(designPar) .* obj.fem.volumes' / obj.options.massLimit;
        end
        
        function gs = constraint2(obj, designPar)
            obj.fem.reassemble(designPar);
            obj.fem.solve();
            gs = obj.fem.displacements' * obj.fem.K_tot * obj.fem.displacements / obj.options.u_max - 1;
        end
        
        function dgsdphi = gradConstraint2(obj, designPar)
            dgsdphi = obj.gradCompliance(designPar) / obj.options.u_max;
        end
    end
    
    methods(Access = protected)
        function dgdphi = gradCompliance(obj, designPar)
            adjointLoads = 2 * obj.fem.K_tot * obj.fem.displacements;
            
            dgdphi = obj.fem.gradChainTerm(adjointLoads);
            
            for e = 1:size(designPar, 2)
                d = designPar(:, e);
                dEdphi = obj.fem.stiffnessDer(d);
                k0 = obj.fem.getElementBaseMatrix(e, 'D');
                u_e = obj.fem.displacements(obj.fem.Edof(:, e), :);
                kT = k0*u_e;
                
                dgdphi(:, e) = dgdphi(:, e) + dEdphi * sum(dot(u_e, kT));
            end
        end
    end
end

