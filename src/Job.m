classdef Job < handle
    properties
        name
        initialSolution
        problemConfiguration
        femConfiguration
        optConf
        result

        problem
    end
    
    methods
        function obj = Job(problem, initial, solverOptions, name)
            if nargin > 3
                obj.name = name;
            end
            obj.problem = problem;
            obj.initialSolution = initial;
            obj.problemConfiguration = problem.getConfiguration();
            obj.femConfiguration = problem.fem.configuration;
            
            obj.optConf.solverOptions = solverOptions;
            obj.optConf.algorithm = solverOptions.algorithm;
        end
        
        function run(obj)
            obj.optConf.solverOptions.min_objective = @(varargin) obj.problem.nlopt_objective(varargin{:});
            obj.optConf.solverOptions.lower_bounds = zeros(size(obj.problem.fem.designPar));
            obj.optConf.solverOptions.upper_bounds = ones(size(obj.problem.fem.designPar));
            obj.optConf.solverOptions.fc = {@(varargin) obj.problem.nlopt_constraint1(varargin{:})};
            
            tic;
            [optDesign, fmin, returnCode] = nlopt_optimize(obj.optConf.solverOptions, obj.initialSolution);
            duration = toc;

            obj.result.finalSolution = optDesign;
            obj.result.fmin = fmin;
            obj.result.returnCode = returnCode;
            obj.result.solverTime = duration;
        end
        
        function save(obj, folderPath)
            % Remove intermediate function which may hold figures
            obj.problem.intermediateFunc = [];
            if nargin < 2
                folderPath = '';
            end
            
            Ex = obj.problem.fem.Ex;
            Ey = obj.problem.fem.Ey;
            Ez = obj.problem.fem.Ey;
               
            saveMatrix = [Ex; Ey; Ez; obj.result.finalSolution'];
            pathPrefix = [folderPath, obj.name];
            jobNameMat = [pathPrefix, '.mat'];
            jobNameTxt = [pathPrefix, '.txt'];
            
            for t = 1:obj.problem.fem.timeSteps
                obj.problem.fem.saveNodeField(pathPrefix + "temperatures" + t, ...
                    obj.problem.fem.temperatures(:, t), "temperatures");
            end
            unfilteredDesign = obj.result.finalSolution;
            filteredDesign = obj.problem.filterParameters(unfilteredDesign);
            obj.problem.fem.saveElementField(pathPrefix + "designFiltered", ...
                filteredDesign, "designFiltered");
            obj.problem.fem.saveElementField(pathPrefix + "designUnfiltered", ...
                unfilteredDesign, "designUnfiltered");
            if isa(obj.femConfiguration.mesh, "StructuredMesh")
                obj.problem.fem.saveElementField(pathPrefix + "designFiltered_POINTS", ...
                    unfilteredDesign, "designUnfiltered", true);
            end
            save(jobNameTxt, 'saveMatrix', '-ascii','-double');
            saveObj = struct(...
                "name", obj.name, ...
                "initialDesign", obj.initialSolution, ...
                "problemConfiguration", obj.problemConfiguration, ...
                "femConfiguration", obj.femConfiguration, ...
                "optimizerConfiguration", obj.optConf, ...
                "result", obj.result ...
                );
            save(jobNameMat, 'saveObj');
        end
        
        function plotResult(obj)
            if obj.problem.fem.spatialDimensions == 3
                warning("The field variable in 3D is not displayable in Matlab");
                return
            end
            ed = obj.problem.fem.getElemTemp(obj.problem.fem.timeSteps-1);
            Ex = obj.problem.fem.Ex;
            Ey = obj.problem.fem.Ey;
               
            figure
            sgtitle(obj.name);
            subplot(1, 2, 1);
            title("Optimal density")
            elfield2(Ex, Ey, obj.result.finalSolution);
            colorbar

            subplot(1, 2, 2);
            title("Terminal temperature distribution");
            elfield2(Ex, Ey, ed);
            colorbar
        end
    end
end