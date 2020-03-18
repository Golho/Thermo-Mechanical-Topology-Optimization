classdef MaxTemperatureProblem < TopOptProblem
    %UNTITLED2 Summary of this class goes here
    %   Detailed explanation goes here

    methods
        function obj = MaxTemperatureProblem(femModel, elementType, options, volumeFraction, intermediateFunc)
            %UNTITLED2 Construct an instance of this class
            %   Detailed explanation goes here
            if nargin ~= 5
                intermediateFunc = [];
            end
            obj = obj@TopOptProblem(femModel, elementType, options, intermediateFunc);
            obj.options.volumeFraction = volumeFraction;
            obj.options.normFactor = 1;
            obj.options.maxSmoothingPar = 15;
        end
        
        function normalize(obj, designPar)
            % Find the maximum temperature for a certain point in design
            % space and use it as a normalization factor
            designPar = reshape(designPar, [], 1);
            if obj.options.filter
                designPar = obj.filterParameters(designPar);
            end
            obj.fem.reassemble(designPar);
            obj.fem.solve();
            obj.options.maxSmoothingPar = 20;
            obj.options.normFactor = max(obj.fem.temperatures, [], 'all');
            obj.options.normFactor = obj.options.normFactor / 100;
        end
        
        function g = objective(obj, designPar)
            designPar = reshape(designPar, [], 1);
            %METHOD1 Summary of this method goes here
            %   Detailed explanation goes here
            if obj.options.filter
                designPar = obj.filterParameters(designPar);
            end
            
            %kappa_2 = obj.options.material_2.kappa;
            obj.fem.reassemble(designPar);
            obj.fem.solve();
            %g = 0;
            %deltaT = obj.fem.tFinal / (obj.fem.timeSteps-1);
            
            % Do smooth maximum over time (aggregate row-wise)
            a = obj.options.maxSmoothingPar;
            g = obj.smax(obj.smax(obj.fem.temperatures(:, 2:end) ...
                / obj.options.normFactor, a, 1), a, 2);
            
            m = max(obj.fem.temperatures(:, 2:end) ...
                / obj.options.normFactor, [], 'all');
            fprintf('Max temp. (norm.): %f\n', m);
            % Detect NaN or Inf from numerical overflow
            if any(isnan(g) | isinf(g), 'all')
                error('The objective function evaluated to either NaN or Inf');
            end
        end
        
        function gs = constraints(obj, designPar)
            designPar = reshape(designPar, [], 1);
            if obj.options.filter
                designPar = obj.filterParameters(designPar);
            end
            designPar = reshape(designPar, [], 1);
            gs(1) = designPar'*obj.fem.mainVolumes / (obj.options.volumeFraction*sum(obj.fem.mainVolumes)) - 1;
        end
        
        function dgdphi = gradObjective(obj, designPar)
            designPar = reshape(designPar, [], 1);
            if obj.options.filter
                designPar = obj.filterParameters(designPar);
            end
            
            %kappa_2 = obj.options.material_2.kappa;
            dgdphi = zeros(length(designPar), 1);
            deltaT = obj.fem.tFinal / (obj.fem.timeSteps-1);
            
            a = obj.options.maxSmoothingPar;

            maxTemp = obj.smax(obj.fem.temperatures(:, 2:end) / obj.options.normFactor, a, 1);
            dWdT = obj.gradSmax(obj.fem.temperatures(:, 2:end) / obj.options.normFactor, a, 1);
            dWdTmax = obj.gradSmax(maxTemp, a, 2);
            adjointLoads = dWdTmax .* dWdT / obj.options.normFactor;
            % Detect NaN or Inf from numerical overflow
            if any(isnan(adjointLoads) | isinf(adjointLoads), 'all')
                error('Some element in the adjoint loads are either NaN or Inf');
            end
            % Adjoint system is solved backward in time
            adjoints = obj.fem.solveAdjoint(adjointLoads);
            
            k0 = obj.fem.getElementBaseMatrix(obj.fem.mainEnod(1, 1), ...
                    obj.elementType, 'D');
            c0 = obj.fem.getElementBaseMatrix(obj.fem.mainEnod(1, 1), ...
                    obj.elementType, 'cp');

            for e = 1:length(designPar)
                edof = obj.fem.mainEnod(e, :);
                T_e = obj.fem.temperatures(edof(2:end), :);
                adjoint_e = adjoints(edof(2:end), :);
                
                
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
            end
            
            if obj.options.filter
                dgdphi = obj.filterGradient(dgdphi);
            end
        end
        
        function dgsdphi = gradConstraints(obj, designPar)
            designPar = reshape(designPar, [], 1);
            dgsdphi(:, 1) = obj.fem.mainVolumes / (obj.options.volumeFraction * sum(obj.fem.mainVolumes));
            
            if obj.options.filter
               dgsdphi(:, 1) = obj.filterGradient(dgsdphi(:, 1)); 
            end
        end
    end
    
    methods(Static)
        function W = smax(X, a, dim)
            W = sum(X.*a.^X, dim) ./ sum(a.^X, dim);
        end
        
        function dW = gradSmax(X, a, dim)
            p = sum(X.*a.^X, dim);
            q = sum(a.^X, dim);
            
            dW = ((1 + X*log(a)).*a.^X) ./ q + ...
                - (a.^X./q).*log(a).*(p./q);
        end
    end
end

