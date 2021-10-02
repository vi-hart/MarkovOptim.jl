module MarkovOptim

using Distributions
using LightGraphs
using LinearAlgebra
using ModelingToolkit

import LightGraphs.SimpleGraphs: SimpleGraphEdge

export RandomParameter
export RateFunction
export MarkovModel, rate_parameters, node_parameters, edge_parameters, rate_function,
    get_graph, steady_states, transition_matrix
export add_edge!, rem_edge!, add_vertex!, rem_vertex!, connect!
export @named

const MTK = ModelingToolkit

include("utils.jl")
include("random_parameter.jl")
include("rate_functions.jl")
include("markov.jl")
include("graphedits.jl")

end
