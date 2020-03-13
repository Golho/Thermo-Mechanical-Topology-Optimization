classdef Job < handle
    properties
        name
        initialSolution
        problemConfiguration
        femConfiguration
        optimizerAlgorithm        
        solverOptions

        problem
        fem
        gmsh
        
        result
        finalSolution
    end
    
    methods
        function obj = Job(problem, initial, solverOptions, name)
            if nargin > 3
                obj.name = name;
            end
            obj.problem = problem;
            obj.initialSolution = initial;
            obj.problemConfiguration = problem.getConfiguration();
            obj.femConfiguration = problem.fem.getConfiguration();
            
            obj.optimizerAlgorithm = solverOptions.algorithm;
            obj.solverOptions = solverOptions;
        end
        
        function run(obj)
            obj.solverOptions.min_objective = @(varargin) obj.problem.nlopt_objective(varargin{:});
            obj.solverOptions.lower_bounds = zeros(size(obj.problem.fem.mainDensities));
            obj.solverOptions.upper_bounds = ones(size(obj.problem.fem.mainDensities));
            obj.solverOptions.fc = {@(varargin) obj.problem.nlopt_constraint1(varargin{:})};
            
            tic;
            [optDesign, fmin, returnCode] = nlopt_optimize(obj.solverOptions, obj.initialSolution);
            duration = toc;

            obj.finalSolution = optDesign;
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
            [Ex, Ey] = obj.problem.fem.getElemTemp(0);
            saveMatrix = [Ex, Ey, obj.finalSolution];
            jobNameMat = [folderPath, obj.name, '.mat'];
            jobNameTxt = [folderPath, obj.name, '.txt'];
            save(jobNameMat, 'obj');
            save(jobNameTxt, 'saveMatrix', '-ascii','-double');
        end
        
        function plotResult(obj)
            [Ex, Ey, ed] = obj.problem.fem.getElemTemp(obj.problem.fem.timeSteps-1);

            figure
            title("Optimal density")
            elfield2(Ex, Ey, obj.finalSolution);
            colorbar

            figure
            title("Terminal temperature distribution");
            elfield2(Ex, Ey, ed);
            colorbar
        end
    end
end