classdef OptimizerMMA < handle
    %OPTIMZERMMA Optimizer of NL optimization problems using the MMA
    %   Detailed explanation goes here
    
    properties
        designPar_old1
        designPar_old2
        designPar
        
        optProblem
        
        terminationFunc
        
        iter = 1;
        upper
        lower
        
        logger
    end
    
    methods
        function obj = OptimizerMMA(optProblem, initialDesignPar, terminationFunc, logger)
            %UNTITLED3 Construct an instance of this class
            %   Detailed explanation goes here
            obj.optProblem = optProblem;
            
            obj.designPar = initialDesignPar;
            obj.upper = ones(size(obj.designPar));
            obj.lower = zeros(size(obj.designPar));
            obj.terminationFunc = terminationFunc;
            
            obj.logger = logger;
        end
        
        function designPar = loop(obj)
            %METHOD1 Summary of this method goes here
            %   Detailed explanation goes here
            while obj.iter == 1 || ~obj.terminationFunc(obj.designPar, obj.designPar_old1)
                obj.logger.info('loop', sprintf('Optimization iteration: %d', obj.iter));
                xmma = obj.updateDesignPar();
                obj.designPar_old2 = obj.designPar_old1;
                obj.designPar_old1 = obj.designPar;
                obj.designPar = xmma;
                obj.iter = obj.iter + 1;
            end
            designPar = obj.designPar;
        end
        
        function xmma = updateDesignPar(obj)
            tic
            g = obj.optProblem.objective(obj.designPar);
            obj.logger.info('updateDesignPar', sprintf('Objective function: \t\t%f', g));
            gs = obj.optProblem.constraints(obj.designPar);
            i = 1;
            for g_i = gs
                obj.logger.info('updateDesignPar', sprintf('Constraint(%d) function: \t%f', i, g_i));
                i = i + 1;
            end
            
            dg = obj.optProblem.gradObjective(obj.designPar);
            dgs = obj.optProblem.gradConstraints(obj.designPar);
            
            m = length(gs);
            n = length(obj.designPar);
            x_min = obj.optProblem.min*ones(size(obj.designPar));
            x_max = obj.optProblem.max*ones(size(obj.designPar));
            
            [xmma, ~, ~, ~, ~, ~, ~, ~, ~, obj.lower, obj.upper] = ...
                mmasub(m, n, obj.iter, obj.designPar, x_min, x_max, ...
                obj.designPar_old1, obj.designPar_old2, g, dg, 0, ...
                gs, dgs', 0, obj.lower,obj.upper, 1, 0, 1000, 0);

            
            endTime = toc;
            obj.logger.info('updateDesignPar', sprintf('Update time: \t\t\t%f seconds', endTime));
        end
    end
end

