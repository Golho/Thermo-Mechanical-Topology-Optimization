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
            [optDesign, fmin, returnCode] = nlopt_optimize(obj.optConf.solverOptions, reshape(obj.initialSolution, 1, []));
            duration = toc;

            obj.result.finalSolution = reshape(optDesign, size(obj.problem.fem.designPar));
            obj.result.finalSolutionFiltered = obj.problem.filterParameters(obj.result.finalSolution);
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
               
            saveMatrix = [Ex; Ey; Ez; obj.result.finalSolution];
            pathPrefix = [folderPath, obj.name];
            jobNameMat = [pathPrefix, '.mat'];
            jobNameTxt = [pathPrefix, '.txt'];
            
            if isa(obj.problem.fem, "HeatFEMBase")
                for t = 1:obj.problem.fem.timeSteps
                    obj.problem.fem.saveNodeField(pathPrefix + "temperatures" + t, ...
                        obj.problem.fem.temperatures(:, t)', "temperatures");
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
            
            if size(obj.result.finalSolution, 1) == 2
                % Switch the place of the first and the second design
                % parameters, to conform with the was Paraview interpret
                % the values
                unfilteredDesign = flipud(unfilteredDesign);
                filteredDesign = flipud(filteredDesign);
            end
            
            optConf_save = obj.optConf;
            % Remove the function handles from being saved
            nosaveFields = ["min_objective", "max_objective", "fc"];
            for field = nosaveFields
                if isfield(optConf_save.solverOptions, field)
                    optConf_save.solverOptions = rmfield(optConf_save.solverOptions, field);
                end
            end

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
                "optimizerConfiguration", optConf_save, ...
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
            if isa(obj.problem.fem, "OptThermoMechStructured")
                subPlots = 3;
            else
                subPlots = 2;
            end
            iSubPlot = 1;
            subplot(subPlots, 1, iSubPlot);
            iSubPlot = iSubPlot + 1;
            title("Optimal density")
            if size(obj.result.finalSolution, 1) == 2
                patchPlot = elfield2(Ex, Ey, filteredDesign(2, :));
                patchPlot.FaceAlpha = "flat";
                patchPlot.FaceVertexAlphaData = filteredDesign(1, :)';
            elseif size(obj.result.finalSolution, 1) == 1
                elfield2(Ex, Ey, filteredDesign);
            end
            colorbar
            caxis([0, 1]);
            
            if allTimeSteps
                timeSteps = 0:obj.problem.fem.timeSteps-1;
            else
                timeSteps = obj.problem.fem.timeSteps-1;
            end
            
            for timeStep = timeSteps
                if isa(obj.problem.fem, "HeatFEMBase")
                    subplot(subPlots, 1, iSubPlot);
                    iSubPlot = iSubPlot + 1;
                    title("Terminal temperature distribution (" + timeStep + ")");
                    ed = obj.problem.fem.getElemTemp(timeStep);
                    elfield2(Ex, Ey, ed);
                    axis equal
                    colorbar
                end
                if isa(obj.problem.fem, "MechFEMBase") || ...
                        isa(obj.problem.fem, "OptThermoMechStructured")
                    subplot(subPlots, 1, iSubPlot);
                    iSubPlot = iSubPlot + 1;
                    title("Deformed geometry (" + timeStep + ")")
                    if isa(obj.problem.fem, "OptThermoMechStructured")
                        ed = obj.problem.fem.mechFEM.getElemDisp(timeStep);
                    else
                        ed = obj.problem.fem.getElemDisp(timeStep);
                    end
                    scaleFactor = 1;
                    newEx = Ex + scaleFactor*ed(1:2:end, :);
                    newEy = Ey + scaleFactor*ed(2:2:end, :);
                    switch size(obj.result.finalSolution, 1)
                        case 1
                            elfield2(newEx, newEy, filteredDesign);
                        case 2
                            patchPlot = elfield2(newEx, newEy, filteredDesign(2, :));
                            patchPlot.FaceAlpha = "flat";
                            patchPlot.FaceVertexAlphaData = filteredDesign(1, :)';
                        otherwise
                            warning("The number of materials can not properly be displayed");
                    end
                    axis equal
                    colorbar
                    caxis([0, 1]);
                end
                drawnow
            end
            
            if ismethod(obj.problem, "plotResults")
                obj.problem.plotResults();
            end
        end
    end
end