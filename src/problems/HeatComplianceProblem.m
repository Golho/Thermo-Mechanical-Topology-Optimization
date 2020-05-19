classdef HeatComplianceProblem < TopOptProblem
    %UNTITLED2 Summary of this class goes here
    %   Detailed explanation goes here


    methods
        function obj = HeatComplianceProblem(femModel, options, volumeFraction, intermediateFunc)
            %UNTITLED2 Construct an instance of this class
            %   Detailed explanation goes here
            if nargin ~= 5
                intermediateFunc = [];
            end
            obj = obj@TopOptProblem(femModel, options, intermediateFunc);
            obj.options.volumeFraction = volumeFraction;
        end
        
        function g = objective(obj, designPar)
            designPar = reshape(designPar, [], 1);
            %METHOD1 Summary of this method goes here
            %   Detailed explanation goes here
            filteredPar = obj.filterParameters(designPar);
            
            obj.fem.reassemble(filteredPar);
            obj.fem.solve();
            
            kappa_2 = obj.options.material_2.Kappa(1);
            deltaT = obj.fem.tFinal / (obj.fem.timeSteps-1);
            KT = deltaT / (obj.fem.tFinal * kappa_2) * ...
                obj.fem.K * obj.fem.temperatures(:, 2:end);
            g = sum(dot(obj.fem.temperatures(:, 2:end), KT));
        end
        
        function dgdphi = gradObjective(obj, designPar)
            designPar = reshape(designPar, [], 1);
            %filteredPar = obj.filterParameters(designPar);
            
            kappa_2 = obj.options.material_2.Kappa(1);
            deltaT = obj.fem.tFinal / (obj.fem.timeSteps-1);

            adjointLoads = 2*deltaT/obj.fem.tFinal * ...
                obj.fem.K*obj.fem.temperatures(:, 2:end) / kappa_2;
            % Adjoint system is solved backward in time
            
            dgdphi = obj.fem.gradChainTerm(adjointLoads);
            
            for e = 1:length(designPar)
                d = designPar(e);
                dkappadphi = obj.fem.conductivityDer(d);
                k0 = obj.fem.getElementBaseMatrix(e, 'D');
                T_e = obj.fem.temperatures(obj.fem.Enod(:, e), :);
                kT = deltaT/obj.fem.tFinal*dkappadphi/kappa_2 * k0*T_e;
                
                dgdphi(e) = dgdphi(e) + sum(dot(T_e, kT));
            end
            
            dgdphi = obj.filterGradient(dgdphi, designPar);
        end
        
        function g = constraint1(obj, designPar)
            designPar = reshape(designPar, [], 1);
            %filteredPar = obj.filterParameters(designPar);

            g = dot(designPar, obj.fem.volumes) / (obj.options.volumeFraction*sum(obj.fem.volumes)) - 1;
        end
        
        function dgdphi = gradConstraint1(obj, designPar)
            designPar = reshape(designPar, [], 1);
            dgdphi = obj.fem.volumes / (obj.options.volumeFraction * sum(obj.fem.volumes));
            
            dgdphi(:) = obj.filterGradient(dgdphi, designPar);
        end
    end
end

