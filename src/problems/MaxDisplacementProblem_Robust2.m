classdef MaxDisplacementProblem_Robust2 < RobustProblem
    %THERMALLYACTUATEDPROBLEM Summary of this class goes here
    %   Detailed explanation goes here

    
    methods
        function obj = MaxDisplacementProblem_Robust2(femModel, options, massLimit, filters, intermediateFunc)
            %THERMALLYACTUATEDPROBLEM Construct an instance of this class
            %   Detailed explanation goes here
            if nargin < 5
                intermediateFunc = [];
            end
            obj = obj@RobustProblem(femModel, options, massLimit, filters, intermediateFunc);
            [obj.options.densityFunc, obj.options.densityDerFunc] = ...
                densitySIMP(obj.options.materials, 1);
        end
        
        function figures = getFigures(obj)
            figures = getFigures@TopOptProblem(obj);
            figures(end + 1) = obj.thresholdFig;
            figures(end + 1) = obj.objectiveFig;
            figures(end + 1) = obj.projectionsFig;
            
        end
        function g = subObjective(obj, designPar)
            I = obj.fem.getDummy("josse");
            
            obj.fem.reassemble(designPar);
            obj.fem.solve();

            g = -sum(dot(I, obj.fem.displacements));
        end
        
        function dgdphi = subGradObjective(obj, designPar)
            I = obj.fem.getDummy("josse");
            
            adjointLoads = -I;
            
            dgdphi = obj.fem.gradChainTerm(adjointLoads);
        end
        
        function g = volumeConstraint(obj, designPar)
            densities = obj.options.densityFunc(designPar);
            g = dot(designPar(1, :), obj.fem.volumes) / obj.dilatedMassLimit - 1;
        end
        
        function dgsdphi = volumeGradConstraint(obj, designPar)
            dgsdphi = obj.fem.volumes' / obj.dilatedMassLimit;
        end
    end
end

