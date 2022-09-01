using Test
using Gurobi
using PiecewiseLinear
prb = create((x->x^2-6x+2), collect(0.0:0.1:6.0), 5, Gurobi.Optimizer);

fit!(prb);

plot(prb);

@test prb.linearization.residual ≈ 0.0998360655737674