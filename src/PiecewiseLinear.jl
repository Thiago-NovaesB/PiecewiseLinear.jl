module PiecewiseLinear

using JuMP
using Gurobi

include("utils.jl")
include("types.jl")
include("model.jl")
include("user_interface.jl")

export create, prepare!, fit!, plot
end # module
