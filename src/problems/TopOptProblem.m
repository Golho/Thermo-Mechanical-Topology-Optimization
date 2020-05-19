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
            'heavisideFilter', [], ...
            'designFilter', [], ...
            'filterRadius', [], ...
            'filterWeightFunction', [], ...
            'material_1', [], ...
            'material_2', [] ...
        );
        
        intermediateFunc
        heaviside
    end
    
    methods(Abstract)
        g = objective(obj, designPar)
        
        dgdphi = gradObjective(obj, designPar)
    end
    
    methods
        function obj = TopOptProblem(femModel, options, intermediateFunc)
            obj.fem = femModel;

            if nargin == 5
                obj.intermediateFunc = intermediateFunc;
            end
            
            if options.heavisideFilter
                obj.heaviside = HeavisideFilter(0.01, 0.5);
            end
            
            % Use indexing of obj.options to force similar structures
            obj.options(1) = options;
            obj.fem.assemble();
        end
        
        function [val, gradient] = nlopt_objective(obj, designPar)
            disp("Calling objective function");
            val = obj.objective(designPar);
            
            if obj.options.heavisideFilter
                obj.heaviside.update(designPar);
            end
            
            if ~isempty(obj.intermediateFunc)
                obj.intermediateFunc(obj.fem, designPar);
            end
            
            if nargout > 1
                gradient = obj.gradObjective(designPar);
            end
        end
        
        function [gs] = nlopt_constraints(obj)
            c = 0;
            gs = {};
            while ismethod(obj, "constraint" + (c+1))
                c = c + 1;
                gs{c} = @(varargin) obj.nlopt_constraint_i(c, varargin{:});
            end
        end
        
        function [val, gradient] = nlopt_constraint_i(obj, iConstraint, designPar)
            val = obj.("constraint" + iConstraint)(designPar);
            fprintf("Calling constraint(%d): %f\n", iConstraint, val);
            if (nargout > 1)
                gradient = obj.("gradConstraint" + iConstraint)(designPar);
            end
        end
        
        function conf = getConfiguration(obj)
            conf = obj.options;
            conf.problemName = class(obj);
        end
        
        function filteredPar = filterParameters(obj, designPar)
            filteredPar = designPar;
            if obj.options.designFilter
                if ~obj.computedWeights
                    obj.weights = obj.fem.computeWeights(...
                        obj.options.filterRadius, ...
                        obj.options.filterWeightFunction);
                    obj.computedWeights = true;
                end
                filteredPar = obj.weights * filteredPar;
            end
            
            if obj.options.heavisideFilter
                filteredPar = obj.heaviside.filter(filteredPar);
            end
        end
        
        function filteredGrad = filterGradient(obj, grad, designPar)
            filteredGrad = grad;
            filteredPar = designPar;
            if obj.options.designFilter
                if ~obj.computedWeights
                    obj.weights = obj.fem.computeWeights(...
                        obj.options.filterRadius, ...
                        obj.options.filterWeightFunction);
                    obj.computedWeights = true;
                end

                filteredGrad = obj.weights * filteredGrad;
                filteredPar = obj.weights * filteredPar;
            end
            
            if obj.options.heavisideFilter
                filteredGrad = obj.heaviside.gradFilter(filteredPar).*filteredGrad;
                filteredPar = obj.heaviside.filter(filteredPar);
            end
            %filteredGrad = 1./max(designPar, 0.01) .* obj.weights * (designPar.*grad);
        end
    end
end

