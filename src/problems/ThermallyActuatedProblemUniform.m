classdef ThermallyActuatedProblemUniform < TopOptProblem
    %THERMALLYACTUATEDPROBLEM Summary of this class goes here
    %   Detailed explanation goes here
    
    methods
        function obj = ThermallyActuatedProblemUniform(femModel, options, volumeFraction, weight, intermediateFunc)
            %THERMALLYACTUATEDPROBLEM Construct an instance of this class
            %   Detailed explanation goes here
            if nargin < 5
                intermediateFunc = [];
            end
            if nargin < 4
                weight = 0;
            end
            obj = obj@TopOptProblem(femModel, options, intermediateFunc);
            obj.options.volumeFraction = volumeFraction;
            obj.options.intermediateWeight = weight;
        end
        
        function g = objective(obj, designPar)
            designPar = reshape(designPar, [], 1);
            %METHOD1 Summary of this method goes here
            %   Detailed explanation goes here
            filteredPar = obj.filterParameters(designPar);
            
            I = obj.fem.getDummy("josse");
            
            obj.fem.reassemble(filteredPar);
            obj.fem.solve();

            g = -sum(dot(I, obj.fem.displacements)) + ...
                obj.options.intermediateWeight*sum(filteredPar.*(1-filteredPar)) / numel(filteredPar);
        end
        
        function dgdphi = gradObjective(obj, designPar)
            designPar = reshape(designPar, [], 1);
            filteredPar = obj.filterParameters(designPar);

            I = obj.fem.getDummy("josse");
            
            adjointLoads_mech = -I;
            
            dgdphi = obj.fem.gradChainTerm(adjointLoads_mech) + ...
                obj.options.intermediateWeight*(1-2*filteredPar) / numel(filteredPar);
            
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

