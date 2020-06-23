classdef HeatComplianceProblem < TopOptProblem
    %UNTITLED2 Summary of this class goes here
    %   Detailed explanation goes here


    methods
        function obj = HeatComplianceProblem(femModel, options, massLimit, intermediateFunc)
            %UNTITLED2 Construct an instance of this class
            %   Detailed explanation goes here
            if nargin ~= 5
                intermediateFunc = [];
            end
            obj = obj@TopOptProblem(femModel, options, intermediateFunc);
            obj.options.massLimit = massLimit;
            [obj.options.densityFunc, obj.options.densityDerFunc] = ...
                densitySIMP(obj.options.materials, ones(size(obj.fem.designPar, 1), 1));
        end
        
        function g = objective(obj, designPar)
            %METHOD1 Summary of this method goes here
            %   Detailed explanation goes here
            
            obj.fem.reassemble(designPar);
            obj.fem.solve();
            
            kappa_2 = obj.options.materials(2).conductivity;
            deltaT = obj.fem.tFinal / (obj.fem.timeSteps-1);

            KT = deltaT / (obj.fem.tFinal * kappa_2) * ...
                obj.fem.K * obj.fem.temperatures(:, 2:end);
            g = sum(dot(obj.fem.temperatures(:, 2:end), KT));
        end
        
        function dgdphi = gradObjective(obj, designPar)
            
            kappa_2 = obj.options.materials(2).conductivity;
            deltaT = obj.fem.tFinal / (obj.fem.timeSteps-1);

            adjointLoads = 2*deltaT / (obj.fem.tFinal * kappa_2) * ...
                obj.fem.K * obj.fem.temperatures(:, 2:end);
            % Adjoint system is solved backward in time
            
            dgdphi = obj.fem.gradChainTerm(adjointLoads);
            
            for e = 1:size(designPar, 2)
                d = designPar(:, e);
                dkappadphi = obj.fem.conductivityDer(d);
                k0 = obj.fem.getElementBaseMatrix(e, 'D');
                T_e = obj.fem.temperatures(obj.fem.Enod(:, e), :);
                kT = deltaT / (obj.fem.tFinal * kappa_2) * k0*T_e;
                
                dgdphi(:, e) = dgdphi(:, e) + dkappadphi * sum(dot(T_e, kT));
            end
        end
        
        function gs = constraint1(obj, designPar)
            densities = obj.options.densityFunc(designPar);
            gs = dot(densities, obj.fem.volumes) / obj.options.massLimit - 1;
        end
        
        function dgsdphi = gradConstraint1(obj, designPar)
            dgsdphi = obj.options.densityDerFunc(designPar) .* obj.fem.volumes' / obj.options.massLimit;
        end
    end
end

