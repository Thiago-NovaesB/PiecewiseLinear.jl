function _create_model_sec!(prb::Problem)

    prb.model = Model(prb.data.solver)
    JuMP.MOI.set(prb.model, JuMP.MOI.Silent(), true)

    _add_variables_sec!(prb)
    _add_expressions_sec!(prb)
    _add_constraints_sec!(prb)
    _objective_function_sec!(prb)
end

function _add_variables_sec!(prb::Problem)
    model = prb.model
    data = prb.data
    cache = prb.cache

    @variable(model, c[1:cache.N,1:Int(cache.N*(cache.N-1)/2)], Bin)
    @variable(model, u[1:Int(cache.N*(cache.N-1)/2)], Bin)
    @variable(model, y[1:cache.N])
    if data.norm == L1
        @variable(model, 0 <= error[1:cache.N])
    end 
    if data.norm == Linf
        @variable(model, max_error)
    end 
end

function _add_expressions_sec!(prb::Problem)
    model = prb.model
    data = prb.data
    cache = prb.cache
    c = model[:c]

    @expression(model, times[n = 1:Int(cache.N*(cache.N-1)/2)], sum(c[:,n]))
    @expression(model, accumulate[n = 1:cache.N], sum(p*c[n,p] for p in 1:Int(cache.N*(cache.N-1)/2)))
end


function _add_constraints_sec!(prb::Problem)
    model = prb.model
    data = prb.data
    cache = prb.cache
    c = model[:c]
    y = model[:y]
    times = model[:times]
    u = model[:u]
    accumulate = model[:accumulate]
    
    @constraint(model, planes[n = 1:cache.N], y[n] == sum((cache.a[i]*data.grid[n]+cache.b[i])*(c[n,i]) for i in 1:Int(cache.N*(cache.N-1)/2)))
    @constraint(model, choose[n = 1:cache.N], sum(c[n,:]) == 1)
    @constraint(model, transition[n = 1:cache.N-1], accumulate[n] <= accumulate[n+1])
    @constraint(model, limit, sum(u) == data.pieces)
    @constraint(model, investigate[n = 1:Int(cache.N*(cache.N-1)/2)], u[n] >= times[n]/Int(cache.N*(cache.N-1)/2))
    if data.norm == L1
        error = model[:error]
        @constraint(model, residual_up[n = 1:cache.N], error[n] >= y[n]-cache.f[n])
        @constraint(model, residual_dw[n = 1:cache.N], error[n] >= cache.f[n]-y[n])
    elseif data.norm == Linf
        max_error = model[:max_error]
        @constraint(model, residual_up[n = 1:cache.N], max_error >= y[n]-cache.f[n])
        @constraint(model, residual_dw[n = 1:cache.N], max_error >= cache.f[n]-y[n])
    end
end

function _objective_function_sec!(prb::Problem)
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
