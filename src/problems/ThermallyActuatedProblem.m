classdef ThermallyActuatedProblem < TopOptProblem
    %THERMALLYACTUATEDPROBLEM Summary of this class goes here
    %   Detailed explanation goes here
    
    methods
        function obj = ThermallyActuatedProblem(femModel, options, volumeFraction, intermediateFunc)
            %THERMALLYACTUATEDPROBLEM Construct an instance of this class
            %   Detailed explanation goes here
            if nargin < 4
                intermediateFunc = [];
            end
            obj = obj@TopOptProblem(femModel, options, intermediateFunc);
            obj.options.volumeFraction = volumeFraction;
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
        
        function gs = constraint1(obj, designPar)
            %designPar = reshape(designPar, [], 1);
            %designPar(:) = obj.filterParameters(designPar);
            gs = dot(designPar, obj.fem.volumes) / ...
                (obj.options.volumeFraction*sum(obj.fem.volumes)) - 1;
        end
        
        function dgsdphi = gradConstraint1(obj, designPar)
            designPar = reshape(designPar, [], 1);
            %filteredPar = obj.filterParameters(designPar);
            dgsdphi = obj.fem.volumes / (obj.options.volumeFraction * sum(obj.fem.volumes));

            dgsdphi(:) = obj.filterGradient(dgsdphi, designPar);
        end
    end
end

