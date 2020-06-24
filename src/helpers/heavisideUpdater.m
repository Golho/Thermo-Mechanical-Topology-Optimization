function updateFunction = heavisideUpdater(growFactor, counterMod)
    updateFunction = @(rho, beta, counter) updateFun(beta, counter, growFactor, counterMod);
end

function beta = updateFun(beta, counter, growFactor, counterMod)
    if counter > 10 && mod(counter, counterMod) == 0
        beta = beta * growFactor;
    end
    beta = min(beta, 10);
end
