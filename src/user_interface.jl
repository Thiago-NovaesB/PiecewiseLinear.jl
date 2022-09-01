function create(f::Union{Function,Vector{Float64}},
    grid::Vector{Float64},
    pieces::Int,
    solver::Union{DataType,Nothing};
    norm::Norm = L1,
    special::Special = Convex,
    weights::Union{Function,Vector{Float64}} = (x->1.0),
    force_y::Dict{Int,Float64} = Dict{Int,Float64}(),
    )
    prb = Problem()
    
    prb.data.f = f
    prb.data.grid = grid
    prb.data.pieces = pieces
    prb.data.solver = solver
    prb.data.norm = norm
    prb.data.special = special
    prb.data.weights = weights
    prb.data.force_y = force_y
    prepare!(prb)
    return prb
end

function prepare!(prb::Problem)
    _check_consistency(prb)
    _load_cache!(prb)
    _create_model!(prb)
    nothing
end

function _check_consistency(prb::Problem)

    data = prb.data
    @assert(data.solver == Gurobi.Optimizer)
    nothing
end

function fit!(prb::Problem)
    optimize!(prb.model)
    linearization = prb.linearization
    model = prb.model

    linearization.x = value.(model[:x]) 
    linearization.y = value.(model[:y])   
    linearization.a = value.(model[:a])   
    linearization.b = value.(model[:b])   

    linearization.residual = objective_value(model) 
    nothing
end

function plot(prb::Problem)

    data = prb.data
    cache = prb.cache
    linearization = prb.linearization

    a = linearization.a 
    b = linearization.b

    plot(data.grid, cache.f, legend = :outertopleft)
    for p in 1:data.pieces
        plot!(data.grid, a[p]*data.grid.+b[p])
    end
    nothing
end
