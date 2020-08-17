classdef MaxDisplacementProblem_Coupled_Robust2 < RobustProblem
    %THERMALLYACTUATEDPROBLEM Summary of this class goes here
    %   Detailed explanation goes here

    
    methods
        function obj = MaxDisplacementProblem_Coupled_Robust2(femModel, options, massLimit1, massLimit2, filters, intermediateFunc)
            %THERMALLYACTUATEDPROBLEM Construct an instance of this class
            %   Detailed explanation goes here
            if nargin < 6
                intermediateFunc = [];
            end
            obj = obj@RobustProblem(femModel, options, massLimit1, filters, intermediateFunc);
            obj.options.massLimit2 = massLimit2;
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
            %METHOD1 Summary of this method goes here
            %   Detailed explanation goes here
            I = obj.fem.mechFEM.getDummy("josse");
            
            obj.fem.reassemble(designPar);
            obj.fem.solve();

            g = -sum(dot(I, obj.fem.mechFEM.displacements));
        end
        
        function dgdphi = subGradObjective(obj, designPar)
            %filteredPar = obj.filterParameters(designPar);

            I = obj.fem.mechFEM.getDummy("josse");
            
            adjointLoads_therm = zeros(size(obj.fem.temperatures));
            adjointLoads_mech = -I;
            
            dgdphi = obj.fem.gradChainTerm(adjointLoads_therm, adjointLoads_mech);
        end
        
        function g = volumeConstraint(obj, designPar)
            densities = obj.options.densityFunc(designPar);
            g = dot(designPar(1, :), obj.fem.volumes) / obj.dilatedMassLimit - 1;
        end
        
        function dgsdphi = volumeGradConstraint(obj, designPar)
            dgsdphi = zeros(size(designPar));
            dgsdphi(1, :) = obj.fem.volumes / obj.dilatedMassLimit;
        end
        
%         function g = constraint2(obj, designPar)
%             g = dot(designPar(1, :) .* designPar(2, :), obj.fem.volumes) / (sum(designPar(1, :) .* obj.fem.volumes') * obj.options.massLimit2) - 1;
%         end
%         
%         function dgdphi = gradConstraint2(obj, designPar)
%             dgdphi = zeros(size(designPar));
%             dgdphi(1, :) = (sum(designPar(1, :) .* obj.fem.volumes') * designPar(2, :) .* obj.fem.volumes' - dot(designPar(1, :) .* designPar(2, :), obj.fem.volumes') * obj.fem.volumes') / (sum(designPar(1, :) .* obj.fem.volumes')^2 * obj.options.massLimit2);
%             dgdphi(2, :) = designPar(1, :) .* obj.fem.volumes' / (sum(designPar(1, :) .* obj.fem.volumes') * obj.options.massLimit2);
%         end
    end
end

