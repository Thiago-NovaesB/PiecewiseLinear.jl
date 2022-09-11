function create(f::Union{Function,Vector{Float64}},
    grid::Vector{Float64},
    pieces::Int,
    solver::Union{DataType,Nothing};
    rules::Rules = free,
    g::Union{Function,Vector{Float64}} = nothing,
    norm::Norm = L2,
    force_prop::ForceProp = none,
    weights::Union{Function,Vector{Float64}} = (x->1.0),
    force_y::Dict{Int,Float64} = Dict{Int,Float64}(),
    )
    prb = Problem();
    
    prb.data.f = f;
    prb.data.g = g;
    prb.data.grid = grid;
    prb.data.pieces = pieces;
    prb.data.solver = solver;
    prb.data.rules = rules;
    prb.data.norm = norm;
    prb.data.force_prop = force_prop;
    prb.data.weights = weights;
    prb.data.force_y = force_y;
    prepare!(prb);
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
    @assert(data.rules == free || !(force_prop == none))
    nothing
end

function _load_cache!(prb::Problem)
    data = prb.data
    cache = prb.cache

    cache.N = length(data.grid)
    if typeof(data.f) == Vector{Float64}
        cache.f = data.f
    else
        cache.f = data.f.(data.grid)
    end

    if typeof(data.g) == Vector{Float64}
        cache.g = data.g
    elseif typeof(data.g) == Function
        cache.g = data.g.(data.grid)
    else
        cache.g = nothing
    end

    if typeof(data.weights) == Vector{Float64}
        cache.w = data.weights
    else
        cache.w = data.weights.(data.grid)
    end

    if data.rule == sec
        cache.a = [(cache.f[n+1] - cache.f[n])/(data.grid[n+1] - data.grid[n]) for n in 1:cache.N-1]
        cache.b = [cache.f[n] - a[n]*data.grid[n] for n in 1:cache.N-1]
    elseif data.rule == grad
        cache.a = cache.g
        cache.b = [cache.f[n] - a[n]*data.grid[n] for n in 1:cache.N]
    end
end

function _create_model!(prb::Problem)

    rule = prb.data.rule

    if rule == free
        _create_model_free!(prb)
    elseif rule == grad
        _create_model_grad!(prb)
    elseif rule == sec
        _create_model_sec!(prb)
    end
    nothing
end

function fit!(prb::Problem)
    optimize!(prb.model)
    linearization = prb.linearization
    model = prb.model
    data = prb.data
    cache = prb.cache

    linearization.x = data.grid 
    linearization.y = value.(model[:y])

       
    if rule == free
        linearization.a = value.(model[:a])   
        linearization.b = value.(model[:b])
    else
        linearization.a = cache.a[value.(model[:c])]   
        linearization.b = cache.b[value.(model[:c])]   
    end
    
    if data.norm == L1
        linearization.residual = objective_value(model)/cache.N
    elseif data.norm == L2
        linearization.residual = sqrt(objective_value(model))/cache.N
    else
        linearization.residual = objective_value(model)
    end
    nothing
end

function plot(prb::Problem)

    data = prb.data
    cache = prb.cache
    linearization = prb.linearization

    a = linearization.a 
    b = linearization.b
    delta_f = maximum(cache.f) - minimum(cache.f)
    ymax =  maximum(cache.f) + 0.1*delta_f
    ymin =  minimum(cache.f) - 0.1*delta_f
    plt = Plots.plot(data.grid, cache.f, lw=5, legend = :outertopleft, title="Error = $(linearization.residual)", label = "f", ylim=(ymin,ymax ))
    for p in 1:data.pieces
        Plots.plot!(plt, data.grid, a[p]*data.grid.+b[p], label="r$(p)", lw=1)
    end
    return plt
end
