classdef ThermallyActuatedProblem < TopOptProblem
    %THERMALLYACTUATEDPROBLEM Summary of this class goes here
    %   Detailed explanation goes here
    
    methods
        function obj = ThermallyActuatedProblem(femModel, options, massLimit, intermediateFunc)
            %THERMALLYACTUATEDPROBLEM Construct an instance of this class
            %   Detailed explanation goes here
            if nargin < 4
                intermediateFunc = [];
            end
            obj = obj@TopOptProblem(femModel, options, intermediateFunc);
            obj.options.massLimit = massLimit;
            [obj.options.densityFunc, obj.options.densityDerFunc] = ...
                densitySIMP(obj.options.materials, 1);
        end
        
        function g = objective(obj, designPar)
            %METHOD1 Summary of this method goes here
            %   Detailed explanation goes here
            I = obj.fem.mechFEM.getDummy("josse");
            
            obj.fem.reassemble(designPar);
            obj.fem.solve();

            g = -sum(dot(I, obj.fem.mechFEM.displacements));
        end
        
        function dgdphi = gradObjective(obj, designPar)
            %filteredPar = obj.filterParameters(designPar);

            I = obj.fem.mechFEM.getDummy("josse");
            
            adjointLoads_therm = zeros(size(obj.fem.temperatures));
            adjointLoads_mech = -I;
            
            dgdphi = obj.fem.gradChainTerm(adjointLoads_therm, adjointLoads_mech);
        end
        
        function gs = constraint1(obj, designPar)
            densities = obj.options.densityFunc(designPar);
            gs = dot(densities, obj.fem.volumes) / obj.options.massLimit - 1;
        end
        
        function dgsdphi = gradConstraint1(obj, designPar)
            dgsdphi = obj.options.densityDerFunc(designPar) .* obj.fem.volumes' / obj.options.massLimit;
        end
    end
end

