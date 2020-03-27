function hQUA_4K(ex, ey, D, t; integration_rule = 2)
    if integration_rule == 1
        g1 = 0
        w1 = 2
        # Gauss points
        gp = [g1, g1]
        # Weights
        w = [w1, w1]

    elseif integration_rule == 2
        g1 = 0.577350269189626
        g2 = 0
        w1 = 1
        w2 = 0.888888888888888
        # Gauss points
        gp = [-g1   -g1;
              g1    -g1;
              -g1   g1;
              g1    g1]

        # Weights
        w = w1*ones(size(gp))

    else
        @error "The integration rule ($integration_rule) is not yet implemented"
    end
    @inbounds begin
        wp = w[:, 1] .* w[:, 2]

        xi = gp[:, 1]
        eta = gp[:, 2]

        nbr_gp = integration_rule^2
        r2 = 2 * nbr_gp

        dNr = zeros(r2+1, 4)
        dNr[1:2:r2, 1] = -(1 .- eta)/4
        dNr[1:2:r2, 3] = (1 .+ eta)/4
        dNr[2:2:r2, 1] = -(1 .- xi)/4
        dNr[2:2:r2+1, 3] = (1 .+ xi)/4

        dNr[1:2:r2, 2] = (1 .- eta)/4
        dNr[1:2:r2, 4] = -(1 .+ eta)/4
        dNr[2:2:r2, 2] = -(1 .+ xi)/4
        dNr[2:2:r2+1, 4] = (1 .- xi)/4

        Ke = zeros(4, 4)
        B = zeros(MMatrix{2, 4, Float64})
        JT = dNr * [ex ey]

        for i = 1:nbr_gp
            index = [2*i-1, 2*i]
            detJ = abs(det(JT[index, :]))
            B = inv(JT[index, :]) * dNr[index, :]
            Ke += B' * D * B * detJ * wp[i]
        end
    end

    return Ke*t
end

function hQUA_4K(ex::Vector{SVector{4,Float64}}, ey::Vector{SVector{4,Float64}}, D::Matrix{SMatrix{2, 2,Float64}}, t::Float64; integration_rule = 2)
    ex = reshape(ex, (4, 1))
    ey = reshape(ey, (4, 1))
    if integration_rule == 1
        g1 = 0
        w1 = 2
        # Gauss points
        gp = [g1, g1]
        # Weights
        w = [w1, w1]

    elseif integration_rule == 2
        g1 = 1 / sqrt(3)
        w1 = 1
        # Gauss points
        gp = [-g1   -g1;
              g1    -g1;
              -g1   g1;
              g1    g1]

        # Weights
        w = w1*ones(size(gp))

    else
        @error "The integration rule ($integration_rule) is not yet implemented"
    end

    wp = w[:, 1] .* w[:, 2]

    xi = gp[:, 1]
    eta = gp[:, 2]

    nbr_gp = integration_rule^2
    r2 = 2 * nbr_gp

    dNr = zeros(r2+1, 4)
    dNr[1:2:r2, 1] = -(1 .- eta)/4
    dNr[1:2:r2, 3] = (1 .+ eta)/4
    dNr[2:2:r2, 1] = -(1 .- xi)/4
    dNr[2:2:r2+1, 3] = (1 .+ xi)/4

    dNr[1:2:r2, 2] = (1 .- eta)/4
    dNr[1:2:r2, 4] = -(1 .+ eta)/4
    dNr[2:2:r2, 2] = -(1 .+ xi)/4
    dNr[2:2:r2+1, 4] = (1 .- xi)/4

    Ke = zeros(4, 4)
    JT = dNr * [ex ey]

    for i = 1:nbr_gp
        index = [2*i-1, 2*i]
        detJ = abs(det(JT[index, :]))
        B = inv(JT[index, :]) * dNr[index, :]
        Ke += B' * D * B * detJ * wp[i]
    end

    return Ke*t
end

function hQUA_4M(ex::Vector{SVector{4,Float64}}, ey::Vector{SVector{4,Float64}}, rho::Float64, t::Float64)
    area::Float64 = abs(ex[3] - ex[1])*abs(ey[3] - ey[1])

    Me::Matrix{SMatrix{4, 4, Float64}} = [4 2 1 2;
         2 4 2 1;
         1 2 4 2;
         2 1 2 4] * 1/36
    rmul!(Me, t * rho * area)
    return Me
end

function hLIN_2f(ex, ey, load_magnitude, t)
    ex = reshape(ex, 2, 1)
    ey = reshape(ey, 2, 1)

    L = sqrt((ex[2]-ex[1])^2 + (ey[2]-ey[1])^2)

    I = 1/2 * [1; 1];
    fe = t * L * load_magnitude * I;

    return fe
end
