classdef HeavisideFilter < handle
    %HEAVISIDEFILTER Summary of this class goes here
    %   Detailed explanation goes here
    
    properties
        beta_0 % Original beta value (used when saving settings)
        beta % non-linearity parameter (steepness)
        omega % threshold value
        updateFun % Function for updating beta
        
        counter = 0;
    end
    
    methods
        function obj = HeavisideFilter(beta_0, omega, updateFun)
            %HEAVISIDEFILTER Construct an instance of this class
            %   Detailed explanation goes here
            obj.beta_0 = beta_0;
            obj.beta = beta_0;
            obj.omega = omega;
            if nargin == 3
                obj.updateFun = updateFun;
            end
        end
        
        function runUpdate(obj, rho)
            if ~isempty(obj.updateFun)
                obj.beta = obj.updateFun(rho, obj.beta, obj.counter);
            else
                obj.beta = obj.update();
            end
            obj.counter = obj.counter + 1;
        end
       
        
        function filteredRho = filter(obj, rho)
            filteredRho = (tanh(obj.beta*obj.omega) + ...
                tanh(obj.beta*(rho - obj.omega))) ...
            / (tanh(obj.beta*obj.omega) + tanh(obj.beta*(1 - obj.omega)));
        end
        
        function gradient = gradFilter(obj, rho)
            gradient = obj.beta*(1 - tanh(obj.beta*(rho-obj.omega)).^2) / ...
            (tanh(obj.beta*obj.omega) + tanh(obj.beta*(1 - obj.omega)));
        end
    end
    
    methods(Access = protected)
         function beta = update(rho, beta, counter)
            if counter > 10 && mod(counter, 20) == 0
                beta = beta * 1.05;
                if beta < 1
                    beta = beta * 1.5;
                elseif beta > 2
                    beta = beta * 2;
                else
                    beta = beta * 1.05;
                end
            end
            beta = min(beta, 10);
        end
    end
end

