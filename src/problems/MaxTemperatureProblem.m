classdef MaxTemperatureProblem < TopOptProblem
    %UNTITLED2 Summary of this class goes here
    %   Detailed explanation goes here

    methods
        function obj = MaxTemperatureProblem(femModel, options, volumeFraction, intermediateFunc)
            %UNTITLED2 Construct an instance of this class
            %   Detailed explanation goes here
            if nargin ~= 5
                intermediateFunc = [];
            end
            obj = obj@TopOptProblem(femModel, options, intermediateFunc);
            obj.options.volumeFraction = volumeFraction;
            obj.options.normFactor = 1;
            obj.options.maxSmoothingPar = 15;
        end
        
        function normalize(obj, designPar)
            % Find the maximum temperature for a certain point in design
            % space and use it as a normalization factor
            designPar = reshape(designPar, [], 1);
            filteredPar = obj.filterParameters(designPar);
            
            obj.fem.reassemble(filteredPar);
            obj.fem.solve();
            obj.options.maxSmoothingPar = 20;
            obj.options.normFactor = max(obj.fem.temperatures, [], 'all');
            obj.options.normFactor = obj.options.normFactor / 100;
        end
        
        function g = objective(obj, designPar)
            designPar = reshape(designPar, [], 1);
            %METHOD1 Summary of this method goes here
            %   Detailed explanation goes here
            filteredPar = obj.filterParameters(designPar);
            
            %kappa_2 = obj.options.material_2.kappa;
            obj.fem.reassemble(filteredPar);
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
        
        function g = constraint1(obj, designPar)
            designPar = reshape(designPar, [], 1);
            filteredPar = obj.filterParameters(designPar);

            g = filteredPar'*obj.fem.volumes / (obj.options.volumeFraction*sum(obj.fem.volumes)) - 1;
        end
        
        function dgdphi = gradObjective(obj, designPar)
            designPar = reshape(designPar, [], 1);
            %filteredPar = obj.filterParameters(designPar);

            a = obj.options.maxSmoothingPar;

            maxTemp = obj.smax(obj.fem.temperatures(:, 2:end) / obj.options.normFactor, a, 1);
            dWdT = obj.gradSmax(obj.fem.temperatures(:, 2:end) / obj.options.normFactor, a, 1);
            dWdTmax = obj.gradSmax(maxTemp, a, 2);
            adjointLoads = dWdTmax .* dWdT / obj.options.normFactor;
            % Detect NaN or Inf from numerical overflow
            if any(isnan(adjointLoads) | isinf(adjointLoads), 'all')
                error('Some element in the adjoint loads are either NaN or Inf');
            end
            
            dgdphi = obj.fem.gradChainTerm(adjointLoads);
            dgdphi(:) = obj.filterGradient(dgdphi, designPar);
        end
        
        function dgdphi = gradConstraint1(obj, designPar)
            designPar = reshape(designPar, [], 1);
            %filteredPar = obj.filterParameters(designPar);
            
            dgdphi = obj.fem.volumes / (obj.options.volumeFraction * sum(obj.fem.volumes));
            dgdphi(:, 1) = obj.filterGradient(dgdphi, designPar);
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
