import Base: convert, insert!, append!

using SparseArrays

struct GlobalAssembly{T<:Number}
    I   ::  Vector{Int}
    J   ::  Vector{Int}
    V   ::  Vector{T}
    m   ::  Int
    n   ::  Int
end

function GlobalAssembly(m::Int, n::Int)
    return GlobalAssembly{Float64}([], [], [], m, n)
end

function GlobalAssembly(m::Int)
    return GlobalAssembly{Float64}([], [], [], m, 1)
end

function convert(::Type{T} where {T<:AbstractSparseMatrix}, A::GlobalAssembly)
    return sparse(A.I, A.J, A.V, A.m, A.n)
end

function insert!(A::GlobalAssembly{T}, I::Int, J::Int, V::T) where T
    push!(A.I, I)
    push!(A.J, J)
    push!(A.V, V)
end

function insert!(A::GlobalAssembly{T}, Is::Vector{Int}, Js::Vector{Int}, e::AbstractArray{T, 2}) where T
    @assert size(e) == (length(Is), length(Js))
    m, n = size(e)
    for i = 1:m, j = 1:i
        insert!(A, Is[i], Js[j], e[i,j])
        if i != j
            insert!(A, Js[j], Is[i], e[i,j])
        end
    end
    return nothing
end

function insert!(A::GlobalAssembly{T}, Is::Vector{Int}, e::AbstractArray{T, 1}) where T
    @assert size(e) == (length(Is),)
    m, = size(e)
    for i = 1:m
        insert!(A, Is[i], 1, e[i])
    end
    return nothing
end

function insert!(A::GlobalAssembly{T}, Is::Vector{Int}, e::T) where T
    @assert length(Is) == 1
    insert!(A, Is[1], 1, e)
    return nothing
end

function append!(main::GlobalAssembly, appendage::GlobalAssembly)
    @assert (main.m, main.n) == (appendage.m, appendage.n)
    append!(main.I, appendage.I)
    append!(main.J, appendage.J)
    append!(main.V, appendage.V)
end
