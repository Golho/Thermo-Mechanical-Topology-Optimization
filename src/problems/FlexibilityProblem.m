classdef FlexibilityProblem < TopOptProblem
    %FLEXIBILITYPROBLEM
    %   An dummy must be prescribed on the FEM model with the name "output"
    
    methods
        function obj = FlexibilityProblem(femModel, options, volumeFraction, u_max, intermediateFunc)
            %UNTITLED2 Construct an instance of this class
            %   Detailed explanation goes here
            if nargin < 5
                intermediateFunc = [];
            end
            obj = obj@TopOptProblem(femModel, options, intermediateFunc);
            obj.options.volumeFraction = volumeFraction;
            obj.options.u_max = u_max;
        end
        
        function g = objective(obj, designPar)
            designPar = reshape(designPar, [], 1);
            %METHOD1 Summary of this method goes here
            %   Detailed explanation goes here
            filteredPar = obj.filterParameters(designPar);
            
            I = obj.fem.getDummy("josse");
            
            obj.fem.reassemble(filteredPar);
            obj.fem.solve();

            g = - I' * obj.fem.displacements / obj.options.u_max;
        end
        
        function dgdphi = gradObjective(obj, designPar)
            designPar = reshape(designPar, [], 1);
            
            %filteredPar = obj.filterParameters(designPar);

            I = obj.fem.getDummy("josse");
            
            adjointLoads = -I / obj.options.u_max;
            
            dgdphi = obj.fem.gradChainTerm(adjointLoads);
            
            dgdphi(:) = obj.filterGradient(dgdphi, designPar);
        end
        
        function gs = constraint1(obj, designPar)
            designPar = reshape(designPar, [], 1);
            filteredPar = obj.filterParameters(designPar);
            gs = dot(filteredPar, obj.fem.volumes) / (obj.options.volumeFraction*sum(obj.fem.volumes)) - 1;
        end
        
        function gs = constraint2(obj, designPar)
            designPar = reshape(designPar, [], 1);
            filteredPar = obj.filterParameters(designPar);
            obj.fem.reassemble(filteredPar);
            obj.fem.solve();
            gs = obj.fem.displacements' * obj.fem.K_tot * obj.fem.displacements / obj.options.u_max - 1;
        end
        
        function dgsdphi = gradConstraint1(obj, designPar)
            designPar = reshape(designPar, [], 1);
            %filteredPar = obj.filterParameters(designPar);
            dgsdphi = obj.fem.volumes / (obj.options.volumeFraction * sum(obj.fem.volumes));

            dgsdphi(:) = obj.filterGradient(dgsdphi, designPar);
        end
        
        function dgsdphi = gradConstraint2(obj, designPar)
            designPar = reshape(designPar, [], 1);
            filteredPar = obj.filterParameters(designPar);
            
            dgsdphi = obj.gradCompliance(filteredPar) / obj.options.u_max;

            dgsdphi(:) = obj.filterGradient(dgsdphi, designPar);
        end
    end
    
    methods(Access = protected)
        function dgdphi = gradCompliance(obj, designPar)
            designPar = reshape(designPar, [], 1);

            adjointLoads = 2 * obj.fem.K_tot * obj.fem.displacements;
            
            dgdphi = obj.fem.gradChainTerm(adjointLoads);
            
            for e = 1:length(designPar)
                d = designPar(e);
                dEdphi = obj.fem.stiffnessDer(d);
                k0 = obj.fem.getElementBaseMatrix(e, 'D');
                u_e = obj.fem.displacements(obj.fem.Edof(:, e), :);
                kT = dEdphi * k0*u_e;
                
                dgdphi(e) = dgdphi(e) + sum(dot(u_e, kT));
            end
        end
    end
end

