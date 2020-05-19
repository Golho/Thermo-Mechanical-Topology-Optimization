classdef Job < handle
    properties
        name
        initialSolution
        problemConfiguration
        femConfiguration
        optFemConf
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
            obj.optFemConf = problem.fem.optConfiguration;
            
            obj.optConf.solverOptions = solverOptions;
            obj.optConf.algorithm = solverOptions.algorithm;
        end
        
        function run(obj)
            obj.optConf.solverOptions.min_objective = @(varargin) obj.problem.nlopt_objective(varargin{:});
            obj.optConf.solverOptions.lower_bounds = zeros(size(obj.problem.fem.designPar));
            obj.optConf.solverOptions.upper_bounds = ones(size(obj.problem.fem.designPar));
            obj.optConf.solverOptions.fc = obj.problem.nlopt_constraints();
            
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
            
            if isa(obj.problem.fem, "HeatFEMBase")
                for t = 1:obj.problem.fem.timeSteps
                    obj.problem.fem.saveNodeField(pathPrefix + "temperatures" + t, ...
                        obj.problem.fem.temperatures(:, t), "temperatures");
                end
            elseif isa(obj.problem.fem, "MechFEMBase")
                for t = 1:obj.problem.fem.timeSteps
                    displacements = reshape(obj.problem.fem.displacements(:, t), ...
                        obj.problem.fem.spatialDimensions, []);
                    obj.problem.fem.saveNodeVectorField(pathPrefix + "displacements" + t, ...
                        displacements, "displacements");
                end
            else
                warning("The result from the FEM model is not registrered as displayable");
            end
            
            if isa(obj.problem.fem, "OptThermoMechStructured")
                for t = 1:obj.problem.fem.mechFEM.timeSteps
                    displacements = reshape(obj.problem.fem.mechFEM.displacements(:, t), ...
                        obj.problem.fem.mechFEM.spatialDimensions, []);
                    obj.problem.fem.mechFEM.saveNodeVectorField(pathPrefix + "displacements" + t, ...
                        displacements, "displacements");
                end
            end
            
            unfilteredDesign = obj.result.finalSolution;
            filteredDesign = obj.problem.filterParameters(unfilteredDesign);
            obj.problem.fem.saveElementField(pathPrefix + "designFiltered", ...
                filteredDesign, "designFiltered");
            obj.problem.fem.saveElementField(pathPrefix + "designUnfiltered", ...
                unfilteredDesign, "designUnfiltered");
            save(jobNameTxt, 'saveMatrix', '-ascii','-double');
            saveObj = struct(...
                "name", obj.name, ...
                "initialDesign", obj.initialSolution, ...
                "problemConfiguration", obj.problemConfiguration, ...
                "femConfiguration", obj.femConfiguration, ...
                "optFemConfiguration", obj.optFemConf, ...
                "optimizerConfiguration", obj.optConf, ...
                "result", obj.result ...
                );
            save(jobNameMat, 'saveObj');
        end
        
        function plotResult(obj, allTimeSteps)
            if nargin < 2
                allTimeSteps = false;
            end
            if obj.problem.fem.spatialDimensions == 3
                warning("The field variable in 3D is not displayable in Matlab");
                return
            end
            Ex = obj.problem.fem.Ex;
            Ey = obj.problem.fem.Ey;
            filteredDesign = obj.problem.filterParameters(obj.result.finalSolution);
               
            figure
            sgtitle(obj.name);
            subplot(1, 2, 1);
            title("Optimal density")
            elfield2(Ex, Ey, filteredDesign);
            colorbar
            
            if allTimeSteps
                timeSteps = 0:obj.problem.fem.timeSteps-1;
            else
                timeSteps = obj.problem.fem.timeSteps-1;
            end

            subplot(1, 2, 2);
            for timeStep = timeSteps
                if isa(obj.problem.fem, "HeatFEMBase")
                    title("Terminal temperature distribution (" + timeStep + ")");
                    ed = obj.problem.fem.getElemTemp(timeStep);
                    elfield2(Ex, Ey, ed);
                    colorbar
                elseif isa(obj.problem.fem, "MechFEMBase")
                    title("Deformed geometry (" + timeStep + ")")
                    ed = obj.problem.fem.getElemDisp(timeStep);
                    eldisplace2(Ex, Ey, ed, 1);
                    axis equal
                else
                    warning("The result from the FEM model is not registrered as displayable");
                end
                if isa(obj.problem.fem, "OptThermoMechStructured")
                    title("Deformed geometry (" + timeStep + ")")
                    ed = obj.problem.fem.mechFEM.getElemDisp(obj.problem.fem.timeSteps-1);
                    eldisplace2(Ex, Ey, ed, 1);
                    axis equal
                end
                drawnow
            end
            
            if ismethod(obj.problem, "plotResults")
                obj.problem.plotResults();
            end
        end
    end
end