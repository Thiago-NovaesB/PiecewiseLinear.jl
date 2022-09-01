# force_y::Dict{Int,Float64}

function _load_cache!(prb::Problem)
    data = prb.data
    cache = prb.cache

    cache.N = length(data.grid)
    if typeof(data.f) == Vector{Float64}
        cache.f = data.f
    else
        cache.f = data.f.(data.grid)
    end

    if typeof(data.weights) == Vector{Float64}
        cache.w = data.weights
    else
        cache.w = data.weights.(data.grid)
    end
end

function _create_model!(prb::Problem)

    prb.model = Model(prb.data.solver)

    _add_variables!(prb)
    _add_expressions!(prb)
    _add_constraints!(prb)
    _objective_function!(prb)
end

function _add_variables!(prb::Problem)
    model = prb.model
    data = prb.data
    cache = prb.cache

    @variable(model, a[1:data.pieces])
    @variable(model, b[1:data.pieces])
    @variable(model, c[1:cache.N,1:data.pieces], Bin)
    @variable(model, z[1:cache.N])
    if data.norm == L1
        @variable(model, error[1:cache.N])
    end 
end

function _add_expressions!(prb::Problem)
    model = prb.model
    data = prb.data
    cache = prb.cache
    c = model[:c]

    @expression(model, accumulate[n = 1:cache.N], sum(p*c[n,p] for p in 1:data.pieces))
end


function _add_constraints!(prb::Problem)
    model = prb.model
    data = prb.data
    cache = prb.cache
    a = model[:a]
    b = model[:b]
    c = model[:c]
    z = model[:z]
    accumulate = model[:accumulate]
    
    @constraint(model, planes[n = 1:cache.N], z[n] == sum((a[i]*data.grid[n]+b[i])*(c[n,i]) for i in 1:data.pieces))
    @constraint(model, choose[n = 1:cache.N], sum(c[n,:]) == 1)
    @constraint(model, transition[n = 1:cache.N-1], accumulate[n] <= accumulate[n+1])
    if data.norm == L1
        error = model[:error]
        @constraint(model, residual_up[n = 1:cache.N], error[n] >= z[n]-cache.f[n])
        @constraint(model, residual_dw[n = 1:cache.N], error[n] >= cache.f[n]-z[n])
    end
    if data.special == Convex
        @constraint(model, convex[i = 1:data.pieces-1], a[i] <= a[i+1])
    elseif data.special == Concave
        @constraint(model, concave[i = 1:data.pieces-1], a[i] >= a[i+1])
    end
end

function _objective_function!(prb::Problem)
    model = prb.model
    cache = prb.cache
    data = prb.data

    if data.norm == L1
        error = model[:error]
        @objective(model, Min, sum(cache.w .* error))
    else
        z = model[:z]
        @objective(model, Min, sum((z[n]-cache.f[n])^2 for n in 1:cache.N))
    end 
end
