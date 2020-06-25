classdef ThermallyActuatedProblem_Robust < TopOptProblem
    %THERMALLYACTUATEDPROBLEM Summary of this class goes here
    %   Detailed explanation goes here
    
    properties
        maxIndex
        dilatedMassLimit
        
        thresholdFig
        thresholdPlots

        objectiveFig
        objectivePlots
        
        projectionsFig
        projectionsPlots
    end
    
    methods
        function obj = ThermallyActuatedProblem_Robust(femModel, options, massLimit, filters, intermediateFunc)
            %THERMALLYACTUATEDPROBLEM Construct an instance of this class
            %   Detailed explanation goes here
            if nargin < 5
                intermediateFunc = [];
            end
            obj = obj@TopOptProblem(femModel, options, intermediateFunc);
            obj.options.massLimit = massLimit;
            obj.dilatedMassLimit = massLimit;
            obj.options.filters = filters;
            [obj.options.densityFunc, obj.options.densityDerFunc] = ...
                densitySIMP(obj.options.materials, 1);
        end
        
        function figures = getFigures(obj)
            figures = getFigures@TopOptProblem(obj);
            figures(end + 1) = obj.thresholdFig;
            figures(end + 1) = obj.objectiveFig;
            figures(end + 1) = obj.projectionsFig;
            
        end
        
        function g = objective(obj, designPar)
            thresholdDesigns = cell(1, 3);
            gs = zeros(3, 1);
            for i = 1:3
                thresholdDesigns{i} = obj.options.filters(i).filter(designPar);
                gs(i) = obj.subObjective(thresholdDesigns{i});
            end
            
            if obj.options.plot
                if obj.iteration == 1
                    obj.initRobustPlotting(thresholdDesigns, obj.options.filters);
                else
                    if obj.fem.spatialDimensions == 2
                        figure(obj.thresholdFig);
                        sgtitle("Designs at iteration " + obj.iteration);
                        for i = 1:3 
                            plotDesign(obj.fem.Ex, obj.fem.Ey, ...
                                thresholdDesigns{i}, obj.thresholdPlots{i});
                        end
                    end
                    
                    figure(obj.objectiveFig);
                    for i = 1:3 
                        obj.objectivePlots(i).XData(end+1) = obj.iteration;
                        obj.objectivePlots(i).YData(end+1) = gs(i);
                        
                        obj.projectionsPlots{i}.YData = obj.options.filters(i).filter(0:0.01:1);
                    end
                end
            end
            
            % Update the mass limit every 10th iteration, and also the
            % first iteration (thereof the minus one)
            if mod(obj.iteration - 1, 10) == 0
                intermediateDensities = obj.options.densityFunc(thresholdDesigns{2});
                intermediateMass = dot(intermediateDensities, obj.fem.volumes);
                dilatedDensities = obj.options.densityFunc(thresholdDesigns{1});
                dilatedMass = dot(dilatedDensities, obj.fem.volumes);
                obj.dilatedMassLimit = dilatedMass / ...
                    intermediateMass * obj.options.massLimit;
            end
            
            [g, obj.maxIndex] = max(gs);
            fprintf("Max index: %d\n", obj.maxIndex);
        end
        
        function g = subObjective(obj, designPar)
            I = obj.fem.mechFEM.getDummy("josse");
            
            obj.fem.reassemble(designPar);
            obj.fem.solve();

            g = -sum(dot(I, obj.fem.mechFEM.displacements));
        end
        
        function dgdphi = gradObjective(obj, designPar)
            filteredDesignPar = obj.options.filters(obj.maxIndex).filter(designPar);
            dgdphi = obj.subGradObjective(filteredDesignPar) .* ...
                obj.options.filters(obj.maxIndex).gradFilter(designPar);
        end
        
        function dgdphi = subGradObjective(obj, designPar)
            I = obj.fem.mechFEM.getDummy("josse");
            
            adjointLoads_therm = zeros(size(obj.fem.temperatures));
            adjointLoads_mech = -I;
            
            dgdphi = obj.fem.gradChainTerm(adjointLoads_therm, adjointLoads_mech);
        end
        
        function g = constraint1(obj, designPar)
            dilatedDesignPar = obj.options.filters(1).filter(designPar);
            g = obj.subConstraint1(dilatedDesignPar);
        end
        
        function g = subConstraint1(obj, designPar)
            densities = obj.options.densityFunc(designPar);
            g = dot(densities, obj.fem.volumes) / obj.dilatedMassLimit - 1;
        end
        
        function dgsdphi = gradConstraint1(obj, designPar)
            dilatedDesignPar = obj.options.filters(1).filter(designPar);
            dgsdphi = obj.options.densityDerFunc(dilatedDesignPar) .* obj.fem.volumes' / obj.dilatedMassLimit;
            dgsdphi(:) = dgsdphi .* obj.options.filters(1).gradFilter(designPar);
        end
        
        function filteredPar = filterForOutput(obj, designPar)
            % The filtering used before the design is considered finished
            filteredPar = obj.filterParameters(designPar);
            filteredPar = obj.options.filters(2).filter(filteredPar);
        end
    end
    
    methods(Access = protected)
        function initRobustPlotting(obj, thresholdDesigns, filters)
            obj.thresholdFig = figure;
            obj.thresholdPlots = cell(size(thresholdDesigns));
            if obj.fem.spatialDimensions == 2
                i = 1;
                sgtitle("Designs at iteration " + obj.iteration);
                for thresholdDesign = thresholdDesigns
                    subplot(3, 1, i);
                    obj.thresholdPlots{i} = plotDesign(obj.fem.Ex, obj.fem.Ey, thresholdDesign{:});
                    title("Design with projection " + i);
                    i = i + 1;
                end
            end
            
            obj.projectionsFig = figure;
            obj.projectionsPlots = cell(numel(filters), 1);
            for i = 1:numel(filters)
                subplot(numel(filters), 1, i);
                obj.projectionsPlots{i} = plot(0:0.01:1, filters(i).filter(0:0.01:1));
                title("Heaviside projection " + i);
            end
            
            obj.objectiveFig = figure;
            obj.objectivePlots = plot(0, zeros(3, 1));
            title("g_" + 0);
            legend(["Dilated", "Intermediate", "Eroded"]);
        end
    end
end

