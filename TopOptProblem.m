classdef TopOptProblem < handle
    %UNTITLED2 Summary of this class goes here
    %   Detailed explanation goes here
    
    properties
        fem
        elementType
        
        p_kappa = 1;
        p_cp = 1;
        
        kappa_1 = 0.001;
        kappa_2 = 2;
        cp_1 = 0.001;
        cp_2 = 2;
    end
    
    methods
        function obj = TopOptProblem(femModel, elementType)
            %UNTITLED2 Construct an instance of this class
            %   Detailed explanation goes here
            obj.fem = femModel;
            obj.elementType = elementType;
            
            cond = @(phi) obj.kappa_1 + phi^obj.p_kappa*...
                (obj.kappa_2-obj.kappa_1);
            cp = @(phi) obj.cp_1 + phi^obj.p_cp*(obj.cp_2 - obj.cp_1);
            dens = @(phi) 1;
            obj.fem.addInterpFuncs(cond, cp, dens);
            
            obj.fem.assemble();
        end
        
        function g = objective(obj, designPar)
            %METHOD1 Summary of this method goes here
            %   Detailed explanation goes here
            obj.fem.reassemble(designPar);
            obj.fem.solve();
            g = 0;
            deltaT = obj.fem.tFinal / (obj.fem.timeSteps-1);
            for n = (obj.fem.timeSteps-1)
                T_n = obj.fem.temperatures(:, n+1);
                g = g + 1/obj.fem.tFinal * deltaT *...
                    T_n'*obj.fem.K*T_n / obj.kappa_2;
            end
        end
        
        function gs = constraints(obj, designPar)
            gs(1) = designPar'*obj.fem.mainVolumes / (0.4*sum(obj.fem.mainVolumes)) - 1;
        end
        
        function dgdphi = derObjective(obj, designPar)
            dgdphi = zeros(length(designPar), 1);
            deltaT = obj.fem.tFinal / (obj.fem.timeSteps-1);
            % Adjoint system is solved backward in time
            adjoints = zeros(obj.fem.nbrDofs, obj.fem.timeSteps-1);
            
            bd = obj.fem.blockedDofs{end};
            freeDofs = ~any(1:obj.fem.nbrDofs == bd);
            dA = decomposition(obj.fem.A(freeDofs, freeDofs));
            B = obj.fem.B;
            T_n = obj.fem.temperatures(:, end);
            rest = 2*deltaT/obj.fem.tFinal * obj.fem.K*T_n / obj.kappa_2;
            for n = (obj.fem.timeSteps-1):-1:1
                Y = rest;
                bd = obj.fem.blockedDofs{n+1};
                freeDofs = ~any(1:obj.fem.nbrDofs == bd);
                partY = Y(freeDofs);
                adjoints(freeDofs, n) = dA \ partY;
                rest = B*adjoints(:, n);
            end
            for e = 1:length(designPar)
                edof = obj.fem.mainEnod(e, :);
                k0 = obj.fem.getElementBaseMatrix(edof(1), ...
                    obj.elementType, 'D');
                c0 = obj.fem.getElementBaseMatrix(edof(1), ...
                    obj.elementType, 'cp');
                dkappadphi = obj.p_kappa*designPar(e)^(obj.p_kappa-1)*...
                    (obj.kappa_2 - obj.kappa_1);
                dcpdphi = obj.p_cp*designPar(e)^(obj.p_cp-1)*...
                    (obj.cp_2 - obj.cp_1);
                for n = 1:(obj.fem.timeSteps-1)
                    T_ne = obj.fem.temperatures(edof(2:end), n+1);
                    T_n_1e = obj.fem.temperatures(edof(2:end), n);
                    adjoint_e = adjoints(edof(2:end), n);
                    dgdphi(e) = dgdphi(e) + ...
                        -adjoint_e' * ( ...
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
                for n = (obj.fem.timeSteps-1)
                    T_ne = obj.fem.temperatures(edof(2:end), n+1);
                    dgdphi(e) = dgdphi(e) + deltaT/obj.fem.tFinal* ...
                        T_ne'*dkappadphi*k0*T_ne / obj.kappa_2;
                end
            end
        end
        
        function dgsdphi = derConstraints(obj, designPar)
            dgsdphi(:, 1) = obj.fem.mainVolumes / (0.4 * sum(obj.fem.mainVolumes));
        end
    end
end

