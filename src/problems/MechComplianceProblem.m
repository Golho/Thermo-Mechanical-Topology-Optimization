classdef MechComplianceProblem < TopOptProblem
    %MECHCOMPLIANCEPROBLEM Summary of this class goes here
    %   Detailed explanation goes here
    
    methods
        function obj = MechComplianceProblem(femModel, options, volumeFraction, intermediateFunc)
            %UNTITLED2 Construct an instance of this class
            %   Detailed explanation goes here
            if nargin < 4
                intermediateFunc = [];
            end
            obj = obj@TopOptProblem(femModel, options, intermediateFunc);
            obj.options.volumeFraction = volumeFraction;
        end
        
        function g = objective(obj, designPar)
            designPar = reshape(designPar, [], 1);
            %METHOD1 Summary of this method goes here
            %   Detailed explanation goes here
            designPar(:) = obj.filterParameters(designPar);
            
            obj.fem.reassemble(designPar);
            obj.fem.solve();
            
            g = obj.fem.displacements' * obj.fem.K_tot * obj.fem.displacements;
        end
        
        function dgdphi = gradObjective(obj, designPar)
            designPar = reshape(designPar, [], 1);
            filteredPar(:) = obj.filterParameters(designPar);

            adjointLoads = 2 * obj.fem.K_tot * obj.fem.displacements;
            
            dgdphi = obj.fem.gradChainTerm(adjointLoads);
            
            for e = 1:length(filteredPar)
                d = filteredPar(e);
                dEdphi = obj.fem.stiffnessDer(d);
                k0 = obj.fem.getElementBaseMatrix(e, 'D');
                u_e = obj.fem.displacements(obj.fem.Edof(:, e), :);
                kT = dEdphi * k0*u_e;
                
                dgdphi(e) = dgdphi(e) + sum(dot(u_e, kT));
            end
            
            dgdphi(:) = obj.filterGradient(dgdphi, designPar);
        end
        
        function gs = constraint1(obj, designPar)
            designPar = reshape(designPar, [], 1);
            designPar = obj.filterParameters(designPar);
            gs = dot(designPar, obj.fem.volumes) / (obj.options.volumeFraction*sum(obj.fem.volumes)) - 1;
        end
        
        function dgsdphi = gradConstraint1(obj, designPar)
            designPar = reshape(designPar, [], 1);
            dgsdphi = obj.fem.volumes / (obj.options.volumeFraction * sum(obj.fem.volumes));
            
            dgsdphi(:) = obj.filterGradient(dgsdphi, designPar);
        end
    end
end

