classdef HeavisideFilter < handle
    %HEAVISIDEFILTER Summary of this class goes here
    %   Detailed explanation goes here
    
    properties
        beta
        omega
    end
    
    properties
        counter = 0;
    end
    
    methods
        function obj = HeavisideFilter(beta_0, omega)
            %HEAVISIDEFILTER Construct an instance of this class
            %   Detailed explanation goes here
            obj.beta = beta_0;
            obj.omega = omega;
        end
        
        function update(obj, rho)
            if obj.counter > 10 && mod(obj.counter, 5) == 0
                if obj.beta < 1
                    obj.beta = obj.beta * 2;
                else
                    obj.beta = obj.beta * 1.05;
                end
            end
            obj.beta = min(obj.beta, 10);
            
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
end

