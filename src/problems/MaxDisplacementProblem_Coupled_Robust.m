classdef MaxDisplacementProblem_Coupled_Robust < RobustProblem
    %THERMALLYACTUATEDPROBLEM Summary of this class goes here
    %   Detailed explanation goes here

    
    methods
        function obj = MaxDisplacementProblem_Coupled_Robust(femModel, options, massLimit, filters, intermediateFunc)
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
            g = dot(densities, obj.fem.volumes) / obj.dilatedMassLimit - 1;
        end
        
        function dgsdphi = volumeGradConstraint(obj, designPar)
            dgsdphi = obj.options.densityDerFunc(designPar) .* obj.fem.volumes' / obj.dilatedMassLimit;
        end      
    end
end

