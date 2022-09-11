function _create_model_sec!(prb::Problem)

    prb.model = BilevelModel(prb.data.solver, mode = BilevelJuMP.SOS1Mode())
    # JuMP.MOI.set(prb.model, JuMP.MOI.Silent(), true)
    _add_variables!(prb)
    _add_constraints!(prb)
    _objective_function!(prb)
end

function _add_variables_sec!(prb::Problem)
    model = prb.model
    data = prb.data
    cache = prb.cache

    @variable(Lower(model), c[1:cache.N-1], Bin)
    @variable(Lower(model), y[1:cache.N])

    if data.norm == L1
        @variable(Upper(model), 0 <= error[1:cache.N])
    end 
    if data.norm == Linf
        @variable(Upper(model), max_error)
    end 
end

function _add_constraints_sec!(prb::Problem)
    model = prb.model
    data = prb.data
    cache = prb.cache
    c = model[:c]
    y = model[:y]

    a = [(cache.f[n+1] - cache.f[n])/(data.grid[n+1] - data.grid[n]) for n in 1:cache.N-1]
    b = [cache.f[n] - a[n]*data.grid[n] for n in 1:cache.N-1]
    if data.force_prop == convex
        @constraint(Lower(model), planes[n = 1:cache.N, i = 1:cache.N-1], y[n] >= (a[i]*data.grid[n]+b[i])*c[i])
    elseif data.force_prop == concave
        @constraint(Lower(model), planes[n = 1:cache.N, i = 1:cache.N-1], y[n] <= (a[i]*data.grid[n]+b[i])*c[i])
    end
    @constraint(Lower(model), choose[n = 1:cache.N-1], sum(c) == data.pieces)
    if data.norm == L1
        error = model[:error]
        @constraint(Upper(model), residual_up[n = 1:cache.N], error[n] >= y[n]-cache.f[n])
        @constraint(Upper(model), residual_dw[n = 1:cache.N], error[n] >= cache.f[n]-y[n])
    elseif data.norm == Linf
        max_error = model[:max_error]
        @constraint(Upper(model), residual_up[n = 1:cache.N], max_error >= y[n]-cache.f[n])
        @constraint(Upper(model), residual_dw[n = 1:cache.N], max_error >= cache.f[n]-y[n])
    end
end

function _objective_function_sec!(prb::Problem)
    model = prb.model
    cache = prb.cache
    data = prb.data

    if data.force_prop == convex
        @objective(Lower(model), Min, sum(y))
    elseif data.force_prop == concave
        @objective(Lower(model), Max, sum(y))
    end
    if data.norm == L1
        error = model[:error]
        @objective(Upper(model), Min, sum(cache.w[n] * error[n] for n in 1:cache.N))
    elseif data.norm == L2
        y = model[:y]
        @objective(Upper(model), Min, sum(cache.w[n] *(y[n]-cache.f[n])^2 for n in 1:cache.N))
    else 
        max_error = model[:max_error]
        @objective(Upper(model), Min, max_error)
    end 
end
