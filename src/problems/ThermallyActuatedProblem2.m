classdef ThermallyActuatedProblem2 < TopOptProblem
    %THERMALLYACTUATEDPROBLEM Summary of this class goes here
    %   Detailed explanation goes here
    
    methods
        function obj = ThermallyActuatedProblem2(femModel, options, intermediateFunc)
            %THERMALLYACTUATEDPROBLEM Construct an instance of this class
            %   Detailed explanation goes here
            if nargin < 3
                intermediateFunc = [];
            end
            obj = obj@TopOptProblem(femModel, options, intermediateFunc);
        end
        
        function g = objective(obj, designPar)
            designPar = reshape(designPar, [], 1);
            %METHOD1 Summary of this method goes here
            %   Detailed explanation goes here
            designPar = obj.filterParameters(designPar);
            
            I = obj.fem.mechFEM.getDummy("josse");
            
            obj.fem.reassemble(designPar);
            obj.fem.solve();

            g = -sum(dot(I, obj.fem.mechFEM.displacements));
        end
        
        function dgdphi = gradObjective(obj, designPar)
            designPar = reshape(designPar, [], 1);
            %filteredPar = obj.filterParameters(designPar);

            I = obj.fem.mechFEM.getDummy("josse");
            
            adjointLoads_therm = zeros(size(obj.fem.temperatures));
            adjointLoads_mech = -I;
            
            dgdphi = obj.fem.gradChainTerm(adjointLoads_therm, adjointLoads_mech);
            
            dgdphi(:) = obj.filterGradient(dgdphi, designPar);
        end
    end
end

