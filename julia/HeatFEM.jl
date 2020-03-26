include("GlobalAssembly.jl")
include("GmshParser.jl")

using SparseArrays
using LinearAlgebra

function assemble(func, element_blocks, node_coordinates, matrix_size)
    nbr_dofs = matrix_size[1]
    g_size = fill(nbr_dofs, length(matrix_size))
    g = GlobalAssembly(g_size...)
    for block in element_blocks
        for element in block.elements
            node_tags = element.node_tags
            ex = node_coordinates[node_tags, 1]
            ey = node_coordinates[node_tags, 2]
            edof = node_to_dof.(node_tags)
            element_matrix = func(ex, ey, block.element_type)

            edofs = fill(edof, length(matrix_size))
            insert!(g, edofs..., element_matrix)
        end
    end
    return g
end

function node_to_dof(node_nbr)
    node_nbr
end

function element_stiffness(ex, ey, D, thickness, element_type::ElementType)
    if element_type == TRI_3
        return hTRI_3K(ex, ey, D, thickness)
    elseif element_type == QUA_4
        return hQUA_4K(ex, ey, D, thickness)
    else
        @error "The element stiffness matrix for element type $element_type is not yet implemented"
    end
end

function element_mass(ex, ey, rho, thickness, element_type::ElementType)
    if element_type == TRI_3
        return hTRI_3M(ex, ey, rho, thickness)
    elseif element_type == QUA_4
        return hQUA_4M(ex, ey, rho, thickness)
    else
        @error "The element mass matrix for element type $element_type is not yet implemented"
    end
end

function element_load(ex, ey, magnitude, thickness, element_type::ElementType)
    if element_type == PNT
        return magnitude*thickness
    elseif element_type == LIN_2
        return hLIN_2f(ex, ey, magnitude, thickness)
    else
        @error "The element load matrix for element type $element_type is not yet implemented"
    end
end

struct BoundaryCondition
    type            :: AbstractString
    physical_name   :: AbstractString
    value           :: Float64
    time_steps      :: Vector{Int}
end

struct BodyCondition
    type            :: AbstractString
    physical_name   :: AbstractString
end

struct Material
    thickness       :: Float64
    density         :: Float64
    heat_capacity   :: Float64
    thermal_cond    :: Matrix{Float64}
end

mutable struct HeatFEM
    gmsh            :: Gmsh
    t_final         :: Float64
    time_steps      :: Int
    nbr_nodes       :: Int
    nbr_dofs        :: Int
    theta           :: Float64
    temperatures    :: Array{Float64, 2}
    K               :: SparseMatrixCSC{Float64, Int}
    C               :: SparseMatrixCSC{Float64, Int}
    F               :: Array{Float64, 2}
    node_coordinates:: Array{Float64, 2}
    blocked_dofs    :: Vector{Vector{Int}}
    body_conditions :: Vector{BodyCondition}
    boundary_conditions :: Vector{BoundaryCondition}
    material        :: Material
end

function HeatFEM(gmsh::Gmsh, t_final::Int, time_steps::Int)
    return HeatFEM(gmsh, convert(Float64, t_final), time_steps)
end

function HeatFEM(gmsh::Gmsh, t_final::Float64, time_steps::Int)
    gmsh = gmsh
    node_coordinates = get_global_coordinates(gmsh)
    nbr_nodes = size(node_coordinates, 1)

    dofs = node_to_dof.(1:nbr_nodes)
    nbr_dofs = length(dofs)

    temperatures = zeros(nbr_dofs, time_steps)

    K = spzeros(nbr_dofs, nbr_dofs)
    C = spzeros(nbr_dofs, nbr_dofs)
    F = zeros(nbr_dofs, time_steps)

    blocked_dofs = fill([], time_steps)
    body_conditions = Vector{BodyCondition}(undef, 0)
    boundary_conditions = Vector{BoundaryCondition}(undef, 0)

    theta = 1

    material = Material(1, 1, 1, Matrix{Float64}(I, 2, 2))

    return HeatFEM(gmsh, t_final, time_steps, nbr_nodes, nbr_dofs, theta,
        temperatures, K, C, F, node_coordinates, blocked_dofs, body_conditions,
        boundary_conditions, material)
end

function get_matrices(obj::HeatFEM)
    dt = obj.t_final / (obj.time_steps - 1)
    A = obj.C + dt*obj.theta*obj.K
    B = obj.C + dt*(obj.theta - 1)*obj.K
    return A, B
end

function solve_transient!(X, A, B, F, bds)
    nbr_dofs, time_steps = size(X)
    bd = bds[1]
    free_dofs = setdiff(1:nbr_dofs, bd)
    dA = cholesky(A[free_dofs, free_dofs])
    for time_step = 1:time_steps-1
        Y = B * X[:, time_step] + F[:, time_step]
        partY = Y[free_dofs] - A[free_dofs, bd]*X[bd, time_step+1]
        X[free_dofs, time_step + 1] = dA \ partY
    end
end

function apply_body_conditions!(obj::HeatFEM)
    D = obj.material.thermal_cond
    s = obj.material.density * obj.material.heat_capacity
    t = obj.material.thickness

    K = GlobalAssembly(obj.nbr_dofs, obj.nbr_dofs)
    C = GlobalAssembly(obj.nbr_dofs, obj.nbr_dofs)
    F_v = spzeros(obj.nbr_dofs, obj.time_steps)
    initial_temp = zeros(obj.nbr_dofs, 1)
    func1(ex, ey, element_type) = element_stiffness(ex, ey, D, t, element_type)
    func2(ex, ey, element_type) = element_mass(ex, ey, s, t, element_type)

    for bc = obj.body_conditions
        element_blocks = get_physical_element_blocks(obj.gmsh, bc.physical_name)
        if bc.type == "main"
            K_i = assemble(func1, element_blocks, obj.node_coordinates,
                (obj.nbr_dofs, obj.nbr_dofs))
            append!(K, K_i)

            C_i = assemble(func2, element_blocks, obj.node_coordinates,
                (obj.nbr_dofs, obj.nbr_dofs))
            append!(C, C_i)

        elseif bc.type == "load"
            magnitude = bc.value
            func3(ex, ey, element_type) = element_load(ex, ey, magnitude, t, element_type)
            F_vi = assemble(func3, element_blocks, obj.node_coordinates,
                (obj.nbr_dofs, ))
            F_vi = convert(SparseMatrixCSC, F_vi)
            F_v[:, bc.time_steps] .+= F_vi

        elseif bc.type == "initial"
            for element_block = element_blocks
                nodes = element_block.elements.node_tags
                dofs = node_to_dof.(nodes)
                initial_temp[dofs] = bc.value
            end
        end
    end
    obj.C = C
    obj.K = K
    obj.temperatures[:, 1] = initial_temp
end

function apply_boundary_conditions!(obj::HeatFEM)
    D = obj.material.thermal_cond
    s = obj.material.density * obj.material.heat_capacity
    t = obj.material.thickness

    F = spzeros(obj.nbr_dofs, obj.time_steps)
    for bc = obj.boundary_conditions
        element_blocks = get_physical_element_blocks(obj.gmsh, bc.physical_name)
        if bc.type == "Neumann"
            magnitude = bc.value
            function func(ex, ey, element_type)
                element_load(ex, ey, magnitude, t, element_type)
            end
            F_i = assemble(func, element_blocks, obj.node_coordinates, (obj.nbr_dofs,))
            F_i = convert(SparseMatrixCSC, F_i)
            F[:, bc.time_steps] .+= F_i
        elseif bc.type == "Robin"
            @error "The robin boundary conditions are not yet implemented"
        elseif bc.type == "Dirichlet"
            nodes = Set{Int}()
            for element_block = element_blocks
                for element in element_block.elements
                    push!(nodes, element.node_tags...)
                end
            end
            nodes_array = collect(nodes)
            for time_step in bc.time_steps
                append!(obj.blocked_dofs[time_step], nodes_array)
            end
            obj.temperatures[nodes_array, bc.time_steps] .= bc.value
        end
    end

    obj.F += F
end

function solve!(obj::HeatFEM)
    dt = obj.t_final / (obj.time_steps - 1)
    weightedLoads = dt*(obj.theta * obj.F[:, 2:end] + (1 - obj.theta)*obj.F[:, 1:end-1])
    A, B = get_matrices(obj)
    solve_transient!(obj.temperatures, A, B, weightedLoads, obj.blocked_dofs)
end

function hQUA_4K(ex, ey, D, t; integration_rule = 2)
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

function hQUA_4K(ex, ey, D, t; integration_rule = 2)
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

function hQUA_4M(ex, ey, rho, t)
    area = abs(ex[3] - ex[1])*abs(ey[3] - ey[1])

    I = [4 2 1 2;
         2 4 2 1;
         1 2 4 2;
         2 1 2 4] * 1/36
    Me = t * rho * area * I
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
