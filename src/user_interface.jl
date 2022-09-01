function create(f::Union{Function,Vector{Float64}},
    grid::Vector{Float64},
    pieces::Int,
    solver::Union{DataType,Nothing};
    norm::Norm = L1,
    force_prop::ForceProp = convex,
    weights::Union{Function,Vector{Float64}} = (x->1.0),
    force_y::Dict{Int,Float64} = Dict{Int,Float64}(),
    )
    prb = Problem();
    
    prb.data.f = f;
    prb.data.grid = grid;
    prb.data.pieces = pieces;
    prb.data.solver = solver;
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
    linearization.a = value.(model[:a])   
    linearization.b = value.(model[:b])   
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
    plt = Plots.plot(data.grid, cache.f, legend = :outertopleft, title="Error = $(linearization.residual)", label = "f", linewdth=1, ylim=(ymin,ymax ))
    for p in 1:data.pieces
        Plots.plot!(plt, data.grid, a[p]*data.grid.+b[p], label="r$(p)", linewdth=5)
    end
    return plt
end
