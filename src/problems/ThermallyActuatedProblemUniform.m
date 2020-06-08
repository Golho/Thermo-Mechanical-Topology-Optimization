classdef ThermallyActuatedProblemUniform < TopOptProblem
    %THERMALLYACTUATEDPROBLEM Summary of this class goes here
    %   Detailed explanation goes here
    
    methods
        function obj = ThermallyActuatedProblemUniform(femModel, options, massLimit, weight, intermediateFunc)
            %THERMALLYACTUATEDPROBLEM Construct an instance of this class
            %   Detailed explanation goes here
            if nargin < 5
                intermediateFunc = [];
            end
            if nargin < 4
                weight = 0;
            end
            obj = obj@TopOptProblem(femModel, options, intermediateFunc);
            obj.options.massLimit = massLimit;
            obj.options.intermediateWeight = weight;
            [obj.options.densityFunc, obj.options.densityDerFunc] = ...
                densitySIMP(obj.options.materials, ones(size(obj.fem.designPar, 1), 1));
        end
        
        function g = objective(obj, designPar)
            %METHOD1 Summary of this method goes here
            %   Detailed explanation goes here
            
            I = obj.fem.getDummy("josse");
            
            obj.fem.reassemble(designPar);
            obj.fem.solve();

            g = -sum(dot(I, obj.fem.displacements)) + ...
                obj.options.intermediateWeight*sum(designPar.*(1-designPar), 'all') / numel(designPar);
        end
        
        function dgdphi = gradObjective(obj, designPar)
            I = obj.fem.getDummy("josse");
            
            adjointLoads_mech = -I;
            
            dgdphi = obj.fem.gradChainTerm(adjointLoads_mech) + ...
                obj.options.intermediateWeight*(1-2*designPar) / numel(designPar);

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

