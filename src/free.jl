function _create_model_free!(prb::Problem)

    prb.model = Model(prb.data.solver)
    JuMP.MOI.set(prb.model, JuMP.MOI.Silent(), true)

    _add_variables_free!(prb)
    _add_expressions_free!(prb)
    _add_constraints_free!(prb)
    _objective_function_free!(prb)
end

function _add_variables_free!(prb::Problem)
    model = prb.model
    data = prb.data
    cache = prb.cache

    @variable(model, a[1:data.pieces])
    @variable(model, b[1:data.pieces])
    @variable(model, c[1:cache.N,1:data.pieces], Bin)
    @variable(model, y[1:cache.N])
    if data.norm == L1
        @variable(model, 0 <= error[1:cache.N])
    end 
    if data.norm == Linf
        @variable(model, max_error)
    end 
end

function _add_expressions_free!(prb::Problem)
    model = prb.model
    data = prb.data
    cache = prb.cache
    c = model[:c]

    @expression(model, accumulate[n = 1:cache.N], sum(p*c[n,p] for p in 1:data.pieces))
end


function _add_constraints_free!(prb::Problem)
    model = prb.model
    data = prb.data
    cache = prb.cache
    a = model[:a]
    b = model[:b]
    c = model[:c]
    y = model[:y]
    accumulate = model[:accumulate]
    
    @constraint(model, planes[n = 1:cache.N], y[n] == sum((a[i]*data.grid[n]+b[i])*(c[n,i]) for i in 1:data.pieces))
    @constraint(model, choose[n = 1:cache.N], sum(c[n,:]) == 1)
    @constraint(model, transition[n = 1:cache.N-1], accumulate[n] <= accumulate[n+1])
    if data.norm == L1
        error = model[:error]
        @constraint(model, residual_up[n = 1:cache.N], error[n] >= y[n]-cache.f[n])
        @constraint(model, residual_dw[n = 1:cache.N], error[n] >= cache.f[n]-y[n])
    elseif data.norm == Linf
        max_error = model[:max_error]
        @constraint(model, residual_up[n = 1:cache.N], max_error >= y[n]-cache.f[n])
        @constraint(model, residual_dw[n = 1:cache.N], max_error >= cache.f[n]-y[n])
    end
    if data.force_prop == convex
        @constraint(model, force_convex[i = 1:data.pieces-1], a[i] <= a[i+1])
        @constraint(model, planes_convex[n = 1:cache.N, i = 1:data.pieces], y[n] >= sum(a[i]*data.grid[n]+b[i]))
    elseif data.force_prop == concave
        @constraint(model, force_concave[i = 1:data.pieces-1], a[i] >= a[i+1])
        @constraint(model, planes_concave[n = 1:cache.N, i = 1:data.pieces], y[n] <= sum(a[i]*data.grid[n]+b[i]))
    end
end

function _objective_function_free!(prb::Problem)
    model = prb.model
    cache = prb.cache
    data = prb.data

    if data.norm == L1
        error = model[:error]
        @objective(model, Min, sum(cache.w[n] * error[n] for n in 1:cache.N))
    elseif data.norm == L2
        y = model[:y]
        @objective(model, Min, sum(cache.w[n] *(y[n]-cache.f[n])^2 for n in 1:cache.N))
    else 
        max_error = model[:max_error]
        @objective(model, Min, max_error)
    end 
end
