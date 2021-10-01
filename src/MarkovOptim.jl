module MarkovOptim

using Distributions
using LightGraphs
using LinearAlgebra
using ModelingToolkit

export RandomParameter
export RateFunction
export MarkovModel, rate_parameters, node_parameters, edge_parameters, rate_function,
    get_graph, steady_states, transition_matrix

const MTK = ModelingToolkit

include("utils.jl")
include("random_parameter.jl")
include("rate_functions.jl")
include("markov.jl")

end
