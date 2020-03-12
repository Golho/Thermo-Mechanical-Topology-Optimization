classdef HeatComplianceProblem < TopOptProblem
    %UNTITLED2 Summary of this class goes here
    %   Detailed explanation goes here

    methods
        function obj = HeatComplianceProblem(femModel, elementType)
            %UNTITLED2 Construct an instance of this class
            %   Detailed explanation goes here
            material_1 = struct('kappa', 0.5, 'cp', 1);
            material_2 = struct('kappa', 20, 'cp', 1);
            
            options = struct(...
                'p_kappa', 3, ...
                'p_cp', 3, ...
                'filter', true, ...
                'filterRadius', 0.02, ...
                'material_1', material_1, ...
                'material_2', material_2 ...
            );
            obj = obj@TopOptProblem(femModel, elementType, options);
        end
        
        function g = objective(obj, designPar)
            designPar = reshape(designPar, [], 1);
            %METHOD1 Summary of this method goes here
            %   Detailed explanation goes here
            if obj.options.filter
                designPar = obj.filterParameters(designPar);
            end
            
            kappa_2 = obj.options.material_2.kappa;
            obj.fem.reassemble(designPar);
            obj.fem.solve();
            g = 0;
            deltaT = obj.fem.tFinal / (obj.fem.timeSteps-1);
            for n = (obj.fem.timeSteps-1):-1:1
                T_n = obj.fem.temperatures(:, n+1);
                g = g + 1/obj.fem.tFinal * deltaT *...
                    T_n'*obj.fem.K*T_n / kappa_2;
            end
        end
        
        function gs = constraints(obj, designPar)
            designPar = reshape(designPar, [], 1);
            if obj.options.filter
                designPar = obj.filterParameters(designPar);
            end
            designPar = reshape(designPar, [], 1);
            gs(1) = designPar'*obj.fem.mainVolumes / (0.4*sum(obj.fem.mainVolumes)) - 1;
        end
        
        function dgdphi = gradObjective(obj, designPar)
            designPar = reshape(designPar, [], 1);
            if obj.options.filter
                designPar = obj.filterParameters(designPar);
            end
            
            kappa_2 = obj.options.material_2.kappa;
            dgdphi = zeros(length(designPar), 1);
            deltaT = obj.fem.tFinal / (obj.fem.timeSteps-1);

            adjointLoads = 2*deltaT/obj.fem.tFinal * ...
                obj.fem.K*obj.fem.temperatures(:, 2:end) / kappa_2;
            % Adjoint system is solved backward in time
            adjoints = obj.fem.solveAdjoint(adjointLoads);

            for e = 1:length(designPar)
                edof = obj.fem.mainEnod(e, :);
                T_e = obj.fem.temperatures(edof(2:end), :);
                adjoint_e = adjoints(edof(2:end), :);
                
                k0 = obj.fem.getElementBaseMatrix(edof(1), ...
                    obj.elementType, 'D');
                c0 = obj.fem.getElementBaseMatrix(edof(1), ...
                    obj.elementType, 'cp');
                dkappadphi = obj.options.p_kappa*designPar(e)^(obj.options.p_kappa-1)*...
                    (obj.options.material_2.kappa - obj.options.material_1.kappa);
                dcpdphi = obj.options.p_cp*designPar(e)^(obj.options.p_cp-1)*...
                    (obj.options.material_2.cp - obj.options.material_1.cp);
                for n = 1:(obj.fem.timeSteps-1)
                    T_ne = T_e(:, n+1);
                    T_n_1e = T_e(:, n);
                    adjoint_ne = adjoint_e(:, n);
                    dgdphi(e) = dgdphi(e) + ...
                        -adjoint_ne' * ( ...
                            ( ...
                                deltaT*obj.fem.theta*dkappadphi*k0 + ...
                                dcpdphi*c0 ...
                            )*T_ne - ...
                            ( ...
                                deltaT*(obj.fem.theta-1)*dcpdphi*k0 + ...
                                dcpdphi*c0 ...
                            )*T_n_1e ...
                        );
                end
                for n = (obj.fem.timeSteps-1):-1:1
                    T_ne = T_e(:, n+1);
                    dgdphi(e) = dgdphi(e) + deltaT/obj.fem.tFinal* ...
                        T_ne'*dkappadphi*k0*T_ne / kappa_2;
                end
            end
            
            if obj.options.filter
                dgdphi = obj.filterGradient(dgdphi);
            end
        end
        
        function dgsdphi = gradConstraints(obj, designPar)
            designPar = reshape(designPar, [], 1);
            dgsdphi(:, 1) = obj.fem.mainVolumes / (0.4 * sum(obj.fem.mainVolumes));
            
            if obj.options.filter
               dgsdphi(:, 1) = obj.filterGradient(dgsdphi(:, 1)); 
            end
        end
    end
end

