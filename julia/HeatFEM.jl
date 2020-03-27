include("GlobalAssembly.jl")
include("GmshParser.jl")
include("ElementRoutines.jl")

using SparseArrays
using LinearAlgebra

function node_to_dof(node_nbr)
    node_nbr
end

# Dispatch on element type

element_stiffness(::Type{Val{QUA_4}}, ex, ey, D, thickness) =
    hQUA_4K(ex, ey, D, thickness)

element_mass(::Type{Val{QUA_4}}, ex, ey, rho, thickness) =
    hQUA_4M(ex, ey, rho, thickness)

element_load(::Type{Val{LIN_2}}, ex, ey, load, thickness) =
    hLIN_2f(ex, ey, load, thickness)

element_load(::Type{Val{PNT}}, ex, ey, load, thickness) =
    load*thickness

struct BoundaryCondition
    type            :: String
    physical_name   :: String
    element_type    :: ElementType
    connectivity    :: Matrix{Int}
    value           :: Float64
    time_steps      :: Vector{Int}
end

function BoundaryCondition(g::Gmsh, type::AbstractString, physical_name::AbstractString,
    value::Float64, time_steps::AbstractArray{Int})
    element_type = get_element_type(g, physical_name)
    connectivity = get_elements(g, physical_name)
    return BoundaryCondition(type, physical_name, element_type, connectivity, value, time_steps)
end

struct BodyCondition
    type            :: String
    physical_name   :: String
    element_type    :: ElementType
    connectivity    :: Matrix{Int}
    value           :: Float64
end

BodyCondition(g::Gmsh, type::AbstractString, physical_name::AbstractString) =
    BodyCondition(g, type, physical_name, 0.0)

function BodyCondition(g::Gmsh, type::AbstractString, physical_name::AbstractString,
    value::Float64)
    element_type = get_element_type(g, physical_name)
    connectivity = get_elements(g, physical_name)
    return BodyCondition(type, physical_name, element_type, connectivity, value)
end


struct Material
    thickness       :: Float64
    density         :: Float64
    heat_capacity   :: Float64
    thermal_cond    :: Matrix{Float64}
end

mutable struct HeatFEM
    node_coordinates:: Matrix{Float64}
    t_final         :: Float64
    time_steps      :: Int
    theta           :: Float64
    nbr_nodes       :: Int
    nbr_dofs        :: Int
    temperatures    :: Matrix{Float64}
    K               :: SparseMatrixCSC{Float64, Int}
    C               :: SparseMatrixCSC{Float64, Int}
    F               :: Matrix{Float64}
    blocked_dofs    :: Vector{Vector{Int}}
    body_conditions :: Vector{BodyCondition}
    boundary_conditions :: Vector{BoundaryCondition}
    material        :: Material
end

function HeatFEM(node_coordinates::Matrix{Float64}, t_final::Int, time_steps::Int)
    return HeatFEM(node_coordinates, convert(Float64, t_final), time_steps)
end

function HeatFEM(node_coordinates::Matrix{Float64}, t_final::Float64, time_steps::Int)
    nbr_nodes = size(node_coordinates, 1)

    dofs = node_to_dof.(1:nbr_nodes)
    nbr_dofs = length(dofs)

    temperatures = zeros(nbr_dofs, time_steps)

    K = spzeros(nbr_dofs, nbr_dofs)
    C = spzeros(nbr_dofs, nbr_dofs)
    F = spzeros(nbr_dofs, time_steps)

    blocked_dofs = fill([], time_steps)
    body_conditions = Vector{BodyCondition}(undef, 0)
    boundary_conditions = Vector{BoundaryCondition}(undef, 0)

    theta = 1
    material = Material(1, 1, 1, Matrix{Float64}(I, 2, 2))

    return HeatFEM(node_coordinates, t_final, time_steps, theta, nbr_nodes, nbr_dofs,
        temperatures, K, C, F, blocked_dofs, body_conditions,
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
    thickness = obj.material.thickness

    K = GlobalAssembly(obj.nbr_dofs, obj.nbr_dofs)
    C = GlobalAssembly(obj.nbr_dofs, obj.nbr_dofs)
    F_v = spzeros(obj.nbr_dofs, obj.time_steps)
    initial_temp = zeros(obj.nbr_dofs, 1)

    for bc = obj.body_conditions
        if bc.type == "main"
            @time begin
            for e = 1:size(bc.connectivity, 1)
                enod = bc.connectivity[e, 2:end]
                edof = node_to_dof.(enod)
                ex = obj.node_coordinates[enod, 1]
                ey = obj.node_coordinates[enod, 2]
                Ke = element_stiffness(Val{bc.element_type}, ex, ey, D, thickness)
                Ce = element_mass(Val{bc.element_type}, ex, ey, s, thickness)

                insert!(K, edof, edof, Ke)
                insert!(C, edof, edof, Ce)
            end
            end

        elseif bc.type == "load"
            load = bc.value
            f = GlobalAssembly(obj.nbr_dofs)
            for e = 1:size(bc.connectivity, 1)
                enod = bc.connectivity[e, 2:end]
                edof = node_to_dof.(enod)
                ex = obj.node_coordinates[enod, 1]
                ey = obj.node_coordinates[enod, 2]
                fe = element_load(Val{bc.element_type}, ex, ey, load, thickness)

                insert!(f, edof, fe)
            end
            F_v[:, bc.time_steps] .+= f

        elseif bc.type == "initial"
            for e = 1:size(bc.connectivity, 1)
                enod = bc.connectivity[e, 2:end]
                edof = node_to_dof.(enod)
                initial_temp[dofs] = bc.value
            end
        end
    end

    obj.C = C
    obj.K = K
    obj.F += F_v
    obj.temperatures[:, 1] = initial_temp
end

function apply_boundary_conditions!(obj::HeatFEM)
    D = obj.material.thermal_cond
    s = obj.material.density * obj.material.heat_capacity
    thickness = obj.material.thickness

    F = spzeros(obj.nbr_dofs, obj.time_steps)
    for bc = obj.boundary_conditions
        if bc.type == "Neumann"
            load = bc.value
            f = GlobalAssembly(obj.nbr_dofs)
            for e = 1:size(bc.connectivity, 1)
                enod = bc.connectivity[e, 2:end]
                edof = node_to_dof.(enod)
                ex = obj.node_coordinates[enod, 1]
                ey = obj.node_coordinates[enod, 2]
                fe = element_load(Val{bc.element_type}, ex, ey, load, thickness)

                insert!(f, edof, fe)
            end
            f = convert(SparseMatrixCSC, f)
            F[:, bc.time_steps] .+= f

        elseif bc.type == "Robin"
            @error "The robin boundary conditions are not yet implemented"

        elseif bc.type == "Dirichlet"
            nodes = unique(bc.connectivity[:, 2:end])
            for time_step in bc.time_steps
                append!(obj.blocked_dofs[time_step], nodes)
            end
            obj.temperatures[nodes, bc.time_steps] .= bc.value
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
