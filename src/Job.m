classdef Job < handle
    %JOB Topology optimization (TO) job 
    %   The TO job handles the optimization of a TO problem and stores the
    %   results. Together with a JobManager, severals Jobs can be run in
    %   the same go, and optimal design from one job can be pipelined to be the
    %   initial design of another job.
    
    properties
        name                    % Name of the job
        initialDesign           % Initial design for optimization of the TO problem
        problemConfiguration    % Configuration of the TO problem
        femConfiguration        % Configuration of the FEM model
        optFemConf              % Configuration of the TO FEM model
        optConf                 % Configuration of the optimizer
        result                  % Struct holding results and auxiliary information

        problem                 % The TO Problem to optimize
        
        linkedJob               % Link the output of linkedJob to be the initial solution 
                                % of this job
    end
    
    methods
        function obj = Job(problem, initial, solverOptions, name)
            if nargin > 3
                obj.name = name;
            end
            obj.problem = problem;
            obj.initialDesign = initial;
            obj.problemConfiguration = problem.getConfiguration();
            obj.femConfiguration = problem.fem.configuration;
            obj.optFemConf = problem.fem.optConfiguration;
            
            obj.optConf.solverOptions = solverOptions;
            obj.optConf.algorithm = solverOptions.algorithm;
        end
        
        function linkJob(obj, parentJob)
            %LINKJOB Link the output of a job to be the initial design of
            %this job
            obj.linkedJob = parentJob;
        end
        
        function run(obj)
            %RUN Send the TO problem to the optimizer
            obj.optConf.solverOptions.min_objective = @(varargin) obj.problem.nlopt_objective(varargin{:});
            obj.optConf.solverOptions.lower_bounds = zeros(size(obj.problem.fem.designPar));
            obj.optConf.solverOptions.upper_bounds = ones(size(obj.problem.fem.designPar));
            obj.optConf.solverOptions.fc = obj.problem.nlopt_constraints();
            
            if ~isempty(obj.linkedJob)
                % Take the output design of the linked job and use as
                % the initial design
                obj.initialDesign = obj.linkedJob.result.finalSolution;
            end

            tic;
            [optDesign, fmin, returnCode] = nlopt_optimize(obj.optConf.solverOptions, reshape(obj.initialDesign, 1, []));
            duration = toc;

            obj.result.finalSolution = reshape(optDesign, size(obj.problem.fem.designPar));
            obj.result.finalSolutionFiltered = obj.problem.filterForOutput(obj.result.finalSolution);
            obj.result.fmin = fmin;
            obj.result.returnCode = returnCode;
            obj.result.solverTime = duration;
        end
        
        function save(obj, folderPath)
            %SAVE Save the results from the optimization to files
            
            if isempty(obj.result)
                warning("No results to save for job " + obj.name)
                return
            end
            
            % Remove intermediate function which may hold figures
            obj.problem.intermediateFunc = [];
            if nargin < 2
                folderPath = "";
            end
            
            Ex = obj.problem.fem.Ex;
            Ey = obj.problem.fem.Ey;
            Ez = obj.problem.fem.Ey;
               
            pathPrefix = folderPath + obj.name;
            jobNameMat = pathPrefix + ".mat";
            jobNameTxt = pathPrefix + ".txt";
            
            % Initiate a VTK Dumper to add data to
            if isa(obj.problem.fem.mesh, "StructuredMesh")
                dumper = VTKDumper(obj.name, obj.problem.fem.mesh);
            else
                dumper = VTKDumper(obj.name, obj.problem.fem.mesh, obj.problem.fem.elementType);
            end
            
            % Add the state variables data to the dumper
            if isa(obj.problem.fem, "HeatFEMBase")
                for t = 1:obj.problem.fem.timeSteps
                    data = struct(...
                        "name", "Temperature", ...
                        "type", "scalars", ...
                        "data", obj.problem.fem.temperatures(:, t), ...
                        "nbrComponents", 1, ...
                        "timeSteps", t ...
                    );
                    dumper.addData(data, "point");
                end
            elseif isa(obj.problem.fem, "MechFEMBase")
                displacements = zeros(3, obj.problem.fem.nbrNodes);
                for t = 1:obj.problem.fem.timeSteps
                    displacements(1:obj.problem.fem.spatialDimensions, :) = ...
                        reshape(obj.problem.fem.displacements(:, t), ...
                        obj.problem.fem.spatialDimensions, []);
                    data = struct(...
                        "name", "Displacement", ...
                        "type", "vectors", ...
                        "data", displacements, ...
                        "nbrComponents", 3, ...
                        "timeSteps", t ...
                    );
                    dumper.addData(data, "point");
                end
            else
                warning("The result from the FEM model is not registrered as displayable");
            end
            
            if isa(obj.problem.fem, "OptThermoMechStructured")
                displacements = zeros(3, obj.problem.fem.mechFEM.nbrNodes);
                for t = 1:obj.problem.fem.mechFEM.timeSteps
                    displacements(1:obj.problem.fem.mechFEM.spatialDimensions, :) = ...
                        reshape(obj.problem.fem.mechFEM.displacements(:, t), ...
                        obj.problem.fem.mechFEM.spatialDimensions, []);
                    data = struct(...
                        "name", "Displacement", ...
                        "type", "vectors", ...
                        "data", displacements, ...
                        "nbrComponents", 3, ...
                        "timeSteps", t ...
                    );
                    dumper.addData(data, "point");
                end
            end
            
            unfilteredDesign = obj.result.finalSolution;
            filteredDesign = obj.result.finalSolutionFiltered;
            
            if size(obj.result.finalSolution, 1) == 2
                % Switch the place of the first and the second design
                % parameters, to conform with the way Paraview interpret
                % the values
                unfilteredDesign = flipud(unfilteredDesign);
                filteredDesign = flipud(filteredDesign);
            end
            
            % Add the optimal design data to the dumper
            data = struct(...
                "name", "FilteredDesign", ...
                "type", "scalars", ...
                "data", filteredDesign, ...
                "nbrComponents", size(obj.result.finalSolution, 1), ...
                "timeSteps", 1:obj.problem.fem.timeSteps ...
                );
            dumper.addData(data, "cell");
            
            data = struct(...
                "name", "UnfilteredDesign", ...
                "type", "scalars", ...
                "data", unfilteredDesign, ...
                "nbrComponents", size(obj.result.finalSolution, 1), ...
                "timeSteps", 1:obj.problem.fem.timeSteps ...
                );
            dumper.addData(data, "cell");
            
            optConf_save = obj.optConf;
            % Remove the function handles from being saved
            nosaveFields = ["min_objective", "max_objective", "fc"];
            for field = nosaveFields
                if isfield(optConf_save.solverOptions, field)
                    optConf_save.solverOptions = rmfield(optConf_save.solverOptions, field);
                end
            end
            % Save a backup .txt file for the filtered design
            saveMatrix = [Ex; Ey; Ez; filteredDesign];
            save(jobNameTxt, 'saveMatrix', '-ascii','-double');
            
            % Save a .mat file with all the configuration data
            saveObj = struct(...
                "name", obj.name, ...
                "initialDesign", obj.initialDesign, ...
                "problemConfiguration", obj.problemConfiguration, ...
                "femConfiguration", obj.femConfiguration, ...
                "optFemConfiguration", obj.optFemConf, ...
                "optimizerConfiguration", optConf_save, ...
                "result", obj.result ...
                );
            save(jobNameMat, 'saveObj');
            
            % Dump the dumper data to .VTK-files
            dumper.dump(pathPrefix);
            
            % Save figures if any
            i = 1;
            for fig = obj.problem.getFigures()
                savefig(fig, pathPrefix + "_fig" + i + ".fig");
                i = i + 1;
            end
        end
        
        function plotResult(obj, allTimeSteps)
            %PLOTRESULT Plot the final design and corresponding state
            %variables
            if isempty(obj.result)
                warning("No results to plot for job " + obj.name);
                return
            end
            
            if nargin < 2
                allTimeSteps = false;
            end

            if obj.problem.fem.spatialDimensions == 3
                warning("The field variable in 3D is not displayable in Matlab");
                return
            end
            Ex = obj.problem.fem.Ex;
            Ey = obj.problem.fem.Ey;
            filteredDesign = obj.result.finalSolutionFiltered;
            axis equal
               
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
            plotDesign(Ex, Ey, filteredDesign);
            
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
                    plotDesign(newEx, newEy, filteredDesign);
                    axis equal
                end
                drawnow
            end
        end
    end
end