classdef OptThermoMechStructured < OptHeatFEMStructured
    % Add HeatFEM superclass to inherit some convenient properties
    %OPTTHERMOMECHSTRUCTURED Summary of this class goes here
    %   Detailed explanation goes here
    
    properties
        mechFEM
    end
    
    methods
        function obj = OptThermoMechStructured(optMechFEM, varargin)
            %OPTTHERMOMECHSTRUCTURED Construct an instance of this class
            %   Detailed explanation goes here
            obj = obj@OptHeatFEMStructured(varargin{:});
            obj.mechFEM = copy(optMechFEM);
        end
        
        function assemble(obj)
            assemble@OptHeatFEMStructured(obj);
            obj.mechFEM.assemble();
        end
        
        function reassemble(obj, designPar)
            reassemble@OptHeatFEMStructured(obj, designPar);
            obj.mechFEM.reassemble(designPar);
        end
        
        function solve(obj)
            solve@OptHeatFEMStructured(obj);
            obj.mechFEM.setTemperatures(obj.temperatures);
            obj.mechFEM.solve();
        end
        
        function chainGrad = gradChainTerm(obj, adjointLoads_temp, adjointLoads_disp, dT0dx)
            startTime = tic;
            if nargin < 4
                dT0dx = sparse(size(obj.temperatures, 1), ...
                    numel(obj.designPar));
            end
            
            deltaT = obj.tFinal / (obj.timeSteps - 1);
            
            adjoints_disp = obj.mechFEM.solveAdjoint(adjointLoads_disp);
            adjointLoads_temp = adjointLoads_temp + obj.mechFEM.K_thermal' * adjoints_disp;
            adjoints_temp = obj.solveAdjoint(adjointLoads_temp);
            
            chainGrad = (-adjoints_temp(:, 1)' * obj.B * dT0dx)';
            chainGrad = reshape(chainGrad, size(obj.designPar));
            
            k_TT0 = obj.getElementBaseMatrix(1, 'D');
            c0 = obj.getElementBaseMatrix(1, 'cp');
            k_uu0 = obj.mechFEM.getElementBaseMatrix(1, 'D');
            k_uT0 = obj.mechFEM.getElementBaseMatrix(1, 'D-alpha');
            
            enod = obj.Enod(:, 1);
            edof = obj.mechFEM.Edof(:, 1);
            T_e = obj.temperatures(enod, :);
            u_e = obj.mechFEM.displacements(edof, :);
            adjoint_temp_e = adjoints_temp(enod, :);
            adjoint_disp_e = adjoints_disp(edof, :);
            
            adjoint_dR_Tdx = zeros(size(obj.designPar, 1), 1);
            adjoint_dR_udx = zeros(size(obj.designPar, 1), 1);
            
            for e = 1:size(chainGrad, 2)
                enod(:) = obj.Enod(:, e);
                edof(:) = obj.mechFEM.Edof(:, e);
                T_e(:) = obj.temperatures(enod, :);
                u_e(:) = obj.mechFEM.displacements(edof, :);
                adjoint_temp_e(:) = adjoints_temp(enod, :);
                adjoint_disp_e(:) = adjoints_disp(edof, :);
            
                d = obj.designPar(:, e);
                dkappadphi = obj.conductivityDer(d);
                dcpdphi = obj.heatCapacityDer(d);
                dEdphi = obj.mechFEM.stiffnessDer(d);
                dalphadphi = obj.mechFEM.thermalExpDer(d);
                
                adjoint_dR_Tdx(:) = dkappadphi * sum(dot(adjoint_temp_e, ...
                        deltaT*obj.theta*k_TT0*T_e(:, 2:end) - deltaT*(obj.theta-1)*k_TT0*T_e(:, 1:end-1)...
                    )) + ...
                    dcpdphi * sum(dot(adjoint_temp_e, c0 * (T_e(:, 2:end) - T_e(:, 1:end-1))));
                
                adjoint_dR_udx(:) = (...
                    dEdphi * sum(dot(adjoint_disp_e, k_uu0 * u_e)) - ...
                    dalphadphi * sum(dot(adjoint_disp_e, k_uT0 * T_e)) ...
                );

                chainGrad(:, e) = chainGrad(:, e) - adjoint_dR_Tdx - adjoint_dR_udx;
            end
            fprintf("Computed gradient term:\t\t%f secs\n", toc(startTime));
        end
    end
    
    methods(Access = protected)
        function conf = getConfiguration(obj)
            conf{1} = getConfiguration@OptHeatFEMStructured(obj);
            conf{2} = obj.mechFEM.configuration;
        end
        
        function optConf = getOptConfiguration(obj)
            optConf{1} = getOptConfiguration@OptHeatFEMStructured(obj);
            optConf{2} = obj.mechFEM.optConfiguration;
        end
    end
end

