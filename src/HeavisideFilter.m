classdef HeavisideFilter < handle
    %HEAVISIDEFILTER Smooth approximation of a Heaviside projection
    %   This class holds the function and derivative for a smooth
    %   approximation of a Heaviside projection. It also allows gradually
    %   updating the steepness of the approximation.
    
    properties
        beta_0          % Original beta value (used when saving settings)
        beta            % non-linearity parameter (steepness)
        eta             % threshold value
        updateFun       % Function for updating beta
        
        counter = 0;    % Internal counter for updating beta
    end
    
    methods
        function obj = HeavisideFilter(beta_0, eta, updateFun)
            obj.beta_0 = beta_0;
            obj.beta = beta_0;
            obj.eta = eta;
            if nargin == 3
                obj.updateFun = updateFun;
            end
        end
        
        function runUpdate(obj, rho)
            %RUNUPDATE Update the steepness of the approximation
            if ~isempty(obj.updateFun)
                obj.beta = obj.updateFun(rho, obj.beta, obj.counter);
            else
                obj.beta = obj.update();
            end
            obj.counter = obj.counter + 1;
        end
       
        
        function filteredRho = filter(obj, rho)
            %FILTER Use the projection to filter the input rho
            filteredRho = (tanh(obj.beta*obj.eta) + ...
                tanh(obj.beta*(rho - obj.eta))) ...
            / (tanh(obj.beta*obj.eta) + tanh(obj.beta*(1 - obj.eta)));
        end
        
        function gradient = gradFilter(obj, rho)
            %GRADFILTER Calculate the first order derivative of the filter
            gradient = obj.beta*(1 - tanh(obj.beta*(rho-obj.eta)).^2) / ...
            (tanh(obj.beta*obj.eta) + tanh(obj.beta*(1 - obj.eta)));
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

