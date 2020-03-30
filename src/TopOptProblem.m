classdef (Abstract) TopOptProblem < handle
    %UNTITLED Summary of this class goes here
    %   Detailed explanation goes here
    
    properties
        fem
        weights
        computedWeights = false;
        
        designMin = 0;
        designMax = 1;
        
        options = struct(...
            'p_kappa', [], ...
            'p_cp', [], ...
            'filter', [], ...
            'filterRadius', [], ...
            'filterWeightFunction', [], ...
            'material_1', [], ...
            'material_2', [] ...
        );
        
        intermediateFunc
    end
    
    methods(Abstract)
        g = objective(obj, designPar)
        gs = constraints(obj, designPar)
        
        dgdphi = gradObjective(obj, designPar)
        dgsdphi = gradConstraints(obj, designPar)
    end
    
    methods
        function obj = TopOptProblem(femModel, options, intermediateFunc)
            obj.fem = femModel;

            if nargin == 5
                obj.intermediateFunc = intermediateFunc;
            end
            
            % Use indexing of obj.options to force similar structures
            obj.options(1) = options;
            
            kappa_1 = obj.options.material_1.Kappa(1);
            kappa_2 = obj.options.material_2.Kappa(1);
            cp_1 = obj.options.material_1.heatCapacity;
            cp_2 = obj.options.material_2.heatCapacity;
            p_kappa = obj.options.p_kappa;
            p_cp = obj.options.p_cp;
            
            cond = @(phi) kappa_1 + phi^p_kappa * (kappa_2 - kappa_1);
            condDer = @(phi) p_kappa*(phi)^(p_kappa-1) * (kappa_2 - kappa_1);
            
            cp = @(phi) cp_1 + phi^p_cp * (cp_2 - cp_1);
            cpDer = @(phi) p_cp*phi^(p_cp-1) * (cp_2 - cp_1);
            
            obj.fem.addInterpFuncs(cond, condDer, cp, cpDer);
            obj.fem.assemble();
        end
        
        function [val, gradient] = nlopt_objective(obj, designPar)
            disp("Calling objective function");
            val = obj.objective(designPar);
            
            if ~isempty(obj.intermediateFunc)
                obj.intermediateFunc(obj.fem, designPar);
            end
            
            if (nargout > 1)
                gradient = obj.gradObjective(designPar);
            end
        end
        
        function [val, gradient] = nlopt_constraint1(obj, designPar)
            gs = obj.constraints(designPar);
            fprintf('Calling constraint(1): %f\n', gs(1));
            val = gs(1);
            if (nargout > 1)
                dgs = obj.gradConstraints(designPar);
                gradient = dgs(:, 1);
            end
        end
        
        % TODO: handle the case when there are more than 1 constraints
        
        function conf = getConfiguration(obj)
            conf = obj.options;
        end
        
        function filteredPar = filterParameters(obj, designPar)
            if ~obj.computedWeights
                obj.weights = obj.fem.computeWeights(...
                    obj.options.filterRadius, ...
                    obj.options.filterWeightFunction);
                obj.computedWeights = true;
            end
            
            filteredPar = obj.weights * designPar;
        end
        
        function filteredGrad = filterGradient(obj, grad, designPar)
            if ~obj.computedWeights
                obj.weights = obj.fem.computeWeights(...
                    obj.options.filterRadius, ...
                    obj.options.filterWeightFunction);
                obj.computedWeights = true;
            end
            
            filteredGrad = obj.weights * grad;
            %filteredGrad = 1./max(designPar, 0.01) .* obj.weights * (designPar.*grad);
        end
    end
end

