module PiecewiseLinear

using JuMP
using Gurobi
using Plots

include("utils.jl")
include("types.jl")
include("free.jl")
include("grad.jl")
include("sec.jl")
include("user_interface.jl")

export create, prepare!, fit!, plot
end # module
