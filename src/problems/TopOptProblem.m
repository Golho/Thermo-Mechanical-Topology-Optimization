classdef (Abstract) TopOptProblem < handle
    %TOPOPTPROBLEM Base class for all TO problems formulated
    %   This abstract class creates an interface towards the nlopt
    %   optimizer and handles filtering and plotting for all topology
    %   optimization problems.
    
    properties
        fem                     % FEM Model instance
        weights                 % Weighting matrix used for filtering
        computedWeights = false;% Boolean for if the weights are computed
        
        designMin = 0;          % Min value for the design parameters
        designMax = 1;          % Max value for the design parameters
        
        options = struct(...    % Topology optimization options
            'heavisideFilter', [], ...
            'designFilter', [], ...
            'filterRadius', [], ...
            'filterWeightFunction', [], ...
            'materials', [], ...
            'plot', [] ...
        );
        
        intermediateFunc        % Function to be called inside objective function. 
        % Takes a FEM model as first input and design parameters as second
        heaviside               % Heaviside filtering object
        iteration = 0;          % Counter of iterations
    end
    
    properties(Access = protected)
        designFig
        curveFig
        heavisideFig
        
        designPlot
        curvePlots
        heavisidePlot
        
        nbrConstraints
    end
    
    methods(Abstract)
        g = objective(obj, designPar)           % Objective function
        dgdphi = gradObjective(obj, designPar)  % Gradient of the objective function
    end
    
    methods
        function obj = TopOptProblem(femModel, options, intermediateFunc)
            obj.fem = copy(femModel);
            obj.nbrConstraints = numel(obj.nlopt_constraints);

            if nargin == 3
                obj.intermediateFunc = intermediateFunc;
            end
            
            if isa(options.heavisideFilter, "HeavisideFilter")
                obj.heaviside = options.heavisideFilter;
            end
            
            % Use indexing of obj.options to force similar structures
            obj.options(1) = options;
            obj.fem.assemble();
        end
        
        function figures = getFigures(obj)
            %GETFIGURES Get all the figures used for plotting
            figures = [];
            if obj.options.plot
                figures(end + 1) = obj.designFig;
                figures(end + 1) = obj.curveFig;
                if ~isempty(obj.heaviside)
                    figures(end + 1) = obj.heavisideFig;
                end
            end
        end
        
        function [val, gradient] = nlopt_objective(obj, designPar)
            %NLOPT_OBJECTIVE Objective function called by the nlopt
            %optimizer
            %   The input from the nl_opt optimizer is a 1D vector with 
            %   every design parameter. The output should be the objective
            %   function value and gradient.
            obj.iteration = obj.iteration + 1;
            
            % Reshape the design parameters to have as many columns as
            % number of elements
            designPar = reshape(designPar, [], size(obj.fem.Enod, 2));
            filteredPar = obj.filterParameters(designPar);

            disp("Calling objective function");
            val = obj.objective(filteredPar);
            
            if obj.options.plot
                if obj.iteration == 1
                    obj.initPlotting();
                end
                
                obj.plotDesign(obj.filterForOutput(designPar));
                obj.plotG_i(0, val);
            end
            
            if ~isempty(obj.heaviside)
                obj.heaviside.runUpdate(filteredPar);
            end
            
            if ~isempty(obj.intermediateFunc)
                obj.intermediateFunc(obj.fem, filteredPar);
            end
            
            if nargout > 1
                gradient = obj.gradObjective(filteredPar);    
                gradient(:) = obj.filterGradient(gradient, designPar);
                % Reshape the gradient to comply with the 1D vector the
                % optimizer uses
                gradient = reshape(gradient, 1, []);
            end
        end
        
        function [gs] = nlopt_constraints(obj)
            %NLOPT_CONSTRAINTS Constraint functions called by the nlopt
            %optimizer.
            %   The output should be a cell array with every constraint
            %   function
            c = 0;
            gs = {};
            while ismethod(obj, "constraint" + (c+1))
                c = c + 1;
                gs{c} = @(varargin) obj.nlopt_constraint_i(c, varargin{:});
            end
        end
        
        function [val, gradient] = nlopt_constraint_i(obj, iConstraint, designPar)
            %NLOPT_CONSTRAINT_I Constraint function i called by the nlopt
            %optimizer
            
            % Reshape the design parameters to have as many columns as
            % number of elements
            designPar = reshape(designPar, [], size(obj.fem.Enod, 2));
            filteredPar = obj.filterParameters(designPar);
            val = obj.("constraint" + iConstraint)(filteredPar);
            
            if obj.options.plot
                obj.plotG_i(iConstraint, val);
            end
            
            fprintf("Calling constraint(%d): %f\n", iConstraint, val);
            if (nargout > 1)
                gradient = obj.("gradConstraint" + iConstraint)(filteredPar);
                gradient(:) = obj.filterGradient(gradient, designPar);
                % Reshape the gradient to comply with the 1D vector the
                % optimizer uses
                gradient = reshape(gradient, 1, []);
            end
        end
        
        function conf = getConfiguration(obj)
            conf = obj.options;
            conf.problemName = class(obj);
        end
        
        function filteredPar = filterForOutput(obj, designPar)
            % The filtering used before the design is considered finished
            filteredPar = obj.filterParameters(designPar);
        end
        
        function relErrors = testGradients(obj, designPar, h)
            %TESTGRADIENTS Calculate the of the gradients
            %   testGradients(obj, designPar, h) Calculate the analytical
            %   and numerical gradients for the objective function and
            %   constriants for designPar and compute the relative error.
            %   The parameter h decides the step to use in the numerical
            %   gradient.
            if nargin < 3
                h = 1e-8;
            end
            
            % Numerical gradient of the objective function
            dgdphi = numGrad(@obj.nlopt_objective, designPar, h);

            % Analytical gradient of the objective function
            [~, der_g] = obj.nlopt_objective(designPar);
            objError = norm(dgdphi - der_g) / norm(dgdphi);
            c = 1;
            constraintFuncs = obj.nlopt_constraints;
            constraintErrors = zeros(size(constraintFuncs));
            for constraintFunc = constraintFuncs
                % Numerical gradient of the constraint
                dgdphi = numGrad(constraintFunc{:}, designPar, h);
                % Analytical gradient of the constraint
                [~, der_g] = constraintFunc{:}(designPar);
                constraintErrors(c) = norm(dgdphi - der_g) / norm(dgdphi);
                c = c + 1;
            end
            relErrors = [objError, constraintErrors];
        end
    end
    
    methods(Access = protected)
        function filteredPar = filterParameters(obj, designPar)
            % Filter the design with a density filter and Heaviside
            % projection
            filteredPar = designPar;
            if obj.options.designFilter
                if ~obj.computedWeights
                    obj.weights = obj.fem.computeWeights(...
                        obj.options.filterRadius, ...
                        obj.options.filterWeightFunction)';
                    obj.computedWeights = true;
                end
                filteredPar = designPar * obj.weights;
            end
            
            if ~isempty(obj.heaviside)
                filteredPar = obj.heaviside.filter(filteredPar);
            end
        end
        
        function filteredGrad = filterGradient(obj, grad, designPar)
            % Add the factor to the gradient caused by the density
            % filtering and Heaviside projection
            filteredGrad = grad;
            filteredPar = designPar;
            if obj.options.designFilter
                if ~obj.computedWeights
                    obj.weights = obj.fem.computeWeights(...
                        obj.options.filterRadius, ...
                        obj.options.filterWeightFunction)';
                    obj.computedWeights = true;
                end

                filteredGrad = filteredGrad * obj.weights;
                filteredPar = filteredPar * obj.weights;
            end
            
            if ~isempty(obj.heaviside)
                filteredGrad = obj.heaviside.gradFilter(filteredPar).*filteredGrad;
                %filteredPar = obj.heaviside.filter(filteredPar);
            end
        end
        
        function initPlotting(obj)
            obj.designFig = figure;
            if obj.fem.spatialDimensions == 2
                title("Design at iteration 0");
                obj.designPlot = plotDesign(obj.fem.Ex, obj.fem.Ey, obj.fem.designPar);
            end
            obj.curveFig = figure;
            obj.curvePlots = cell(1, 1 + obj.nbrConstraints);
            for i = 1:(1+obj.nbrConstraints)
                subplot(1 + obj.nbrConstraints, 1, i);
                obj.curvePlots{i} = plot(0, 0);
                title("g_" + (i-1));
            end
            
            if ~isempty(obj.heaviside)
                obj.heavisideFig = figure;
                obj.heavisidePlot = plot(0:0.01:1, obj.heaviside.filter(0:0.01:1));
                title("Heaviside projection (beta = " + obj.heaviside.beta);
            end
        end
        
        function plotDesign(obj, filteredPar)
            if obj.fem.spatialDimensions == 2
                figure(obj.designFig);
                title("Design at iteration " + obj.iteration);
                plotDesign(obj.fem.Ex, obj.fem.Ey, filteredPar, obj.designPlot);
            end
            
            if ~isempty(obj.heaviside)
                figure(obj.heavisideFig);
                obj.heavisidePlot.YData = obj.heaviside.filter(0:0.01:1);
                title("Heaviside projection (beta = " + obj.heaviside.beta);
            end
        end
        
        function plotG_i(obj, iFunction, value)
            figure(obj.curveFig);
            subplot(1 + numel(obj.nlopt_constraints), 1, 1)
            obj.curvePlots{iFunction + 1}.XData(end+1) = obj.iteration;
            obj.curvePlots{iFunction + 1}.YData(end+1) = value;
        end
    end
end

