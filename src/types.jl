@enum Norm L1=1 L2=2
@enum Special Convex=1 Concave=2 None=3

@kwdef mutable struct Data
    norm::Norm = L1
    special::Special = Convex
    grid::Vector{Float64}
    f::Union{Function,Vector{Float64}} = Float64[]
    pieces::Int
    solver::Union{DataType,Nothing}=nothing
    force_y::Dict{Int,Float64}
    weights::Union{Function,Vector{Float64}} = Float64[]
end

@kwdef mutable struct Cache

    N::Int
    f::Vector{Float64}
    w::Vector{Float64}

end

@kwdef mutable struct Linearization
    residual::Float64
    x::Vector{Float64}
    y::Vector{Float64}
    a::Vector{Float64}
    b::Vector{Float64}
end

@kwdef mutable struct Problem
    data::Data
    cache::Cache
    linearization::Linearization
    model::JuMP.Model
end
