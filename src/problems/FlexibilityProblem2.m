classdef FlexibilityProblem2 < TopOptProblem
    %FLEXIBILITYPROBLEM
    %   An dummy must be prescribed on the FEM model with the name "output"
    
    methods
        function obj = FlexibilityProblem2(femModel, options, volumeFraction, intermediateFunc)
            %UNTITLED2 Construct an instance of this class
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
            filteredPar = obj.filterParameters(designPar);
            
            I = obj.fem.getDummy("josse");
            
            obj.fem.reassemble(filteredPar);
            obj.fem.solve();

            g =  -(I' * obj.fem.displacements) / ...
                (obj.fem.displacements' * obj.fem.K_tot * obj.fem.displacements);
        end
        
        function dgdphi = gradObjective(obj, designPar)
            designPar = reshape(designPar, [], 1);
            %filteredPar = obj.filterParameters(designPar);

            I = obj.fem.getDummy("josse");
            inputComp = obj.fem.displacements' * obj.fem.K_tot * obj.fem.displacements;
            outputComp = I' * obj.fem.displacements;
            
            adjointLoads = (inputComp * ...
                -I + ...
                outputComp * 2 * obj.fem.K_tot * obj.fem.displacements)...
                / inputComp^2;
            
            dgdphi = obj.fem.gradChainTerm(adjointLoads);
            
            for e = 1:length(designPar)
                d = designPar(e);
                dEdphi = obj.fem.stiffnessDer(d);
                k0 = obj.fem.getElementBaseMatrix(e, 'D');
                u_e = obj.fem.displacements(obj.fem.Edof(:, e), :);
                kT = dEdphi * k0*u_e * -outputComp / inputComp^2;
                
                dgdphi(e) = dgdphi(e) + sum(dot(u_e, kT));
            end
            
            
            dgdphi = obj.filterGradient(dgdphi, designPar);
        end
        
        function gs = constraint1(obj, designPar)
            designPar = reshape(designPar, [], 1);
            filteredPar = obj.filterParameters(designPar);
            gs = dot(filteredPar, obj.fem.volumes) / (obj.options.volumeFraction*sum(obj.fem.volumes)) - 1;
        end
        
        function dgsdphi = gradConstraint1(obj, designPar)
            designPar = reshape(designPar, [], 1);
            %filteredPar = obj.filterParameters(designPar);

            dgsdphi = obj.fem.volumes / (obj.options.volumeFraction * sum(obj.fem.volumes));
            dgsdphi(:) = obj.filterGradient(dgsdphi, designPar);
        end
    end
end

