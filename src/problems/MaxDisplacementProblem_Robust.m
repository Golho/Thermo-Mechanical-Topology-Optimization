classdef MaxDisplacementProblem_Robust < RobustProblem
    %THERMALLYACTUATEDPROBLEM Summary of this class goes here
    %   Detailed explanation goes here

    
    methods
        function obj = MaxDisplacementProblem_Robust(femModel, options, massLimit, filters, intermediateFunc)
            %THERMALLYACTUATEDPROBLEM Construct an instance of this class
            %   Detailed explanation goes here
            if nargin < 5
                intermediateFunc = [];
            end
            obj = obj@RobustProblem(femModel, options, massLimit, filters, intermediateFunc);
            [obj.options.densityFunc, obj.options.densityDerFunc] = ...
                densitySIMP(obj.options.materials, 1);
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
            g = dot(densities, obj.fem.volumes) / obj.dilatedMassLimit - 1;
        end
        
        function dgsdphi = volumeGradConstraint(obj, designPar)
            dgsdphi = obj.options.densityDerFunc(designPar) .* obj.fem.volumes' / obj.dilatedMassLimit;
        end
    end
end

