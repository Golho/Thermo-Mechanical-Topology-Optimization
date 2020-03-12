classdef (Abstract) TopOptProblem < handle
    %UNTITLED Summary of this class goes here
    %   Detailed explanation goes here
    
    properties
        fem
        elementType
        weights
        computedWeights = false;
        
        designMin = 0;
        designMax = 1;
        
        options = struct(...
            'p_kappa', [], ...
            'p_cp', [], ...
            'filter', [], ...
            'filterRadius', [], ...
            'material_1', [], ...
            'material_2', [] ...
        );
    end
    
    methods(Abstract)
        g = objective(obj, designPar)
        gs = constraints(obj, designPar)
        
        dgdphi = gradObjective(obj, designPar)
        dgsdphi = gradConstraints(obj, designPar)
    end
    
    methods
        function obj = TopOptProblem(femModel, elementType, options)
            obj.fem = femModel;
            obj.elementType = elementType;
            
            % Explicitly state all the struct properties to ensure the
            % options are valid
            if ~obj.optionsMatching(options)
                error('The passed options does not match the required structure. See the abstract class "TopOptProblem" for the correct structure');
            end
            
            obj.options = options;
            
            kappa_1 = obj.options.material_1.kappa;
            kappa_2 = obj.options.material_2.kappa;
            cp_1 = obj.options.material_1.cp;
            cp_2 = obj.options.material_2.cp;
            p_kappa = obj.options.p_kappa;
            p_cp = obj.options.p_cp;
            
            cond = @(phi) kappa_1 + phi^p_kappa*...
                (kappa_2 - kappa_1);
            cp = @(phi) cp_1 + phi^p_cp*(cp_2 - cp_1);
            dens = @(phi) 1;
            
            obj.fem.addInterpFuncs(cond, cp, dens);
            obj.fem.assemble();
        end
        
        function [val, gradient] = nlopt_objective(obj, designPar)
            disp("Calling objective function");
            val = obj.objective(designPar);
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
        
        function match = optionsMatching(obj, options)
            match = true;
            for fieldname = fieldnames(obj.options)
                if ~isfield(options, fieldname)
                    match = false;
                    break;
                end
            end
        end
        
        function filteredPar = filterParameters(obj, designPar)
            if ~obj.computedWeights
                obj.weights = obj.fem.computeMainWeights(obj.options.filterRadius);
                obj.computedWeights = true;
            end
            
            filteredPar = obj.weights * designPar;
        end
        
        function filteredGrad = filterGradient(obj, grad, designPar)
            if ~obj.computedWeights
                obj.weights = obj.fem.computeMainWeights(obj.options.filterRadius);
                obj.computedWeights = true;
            end
            
            filteredGrad = obj.weights * grad;
            %filteredGrad = 1./max(designPar, 0.01) .* obj.weights * (designPar.*grad);
        end
    end
end

