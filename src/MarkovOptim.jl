module MarkovOptim

using Distributions
using LightGraphs
using LinearAlgebra
using ModelingToolkit

import LightGraphs.SimpleGraphs: SimpleGraphEdge
import ModelingToolkit: AbstractTimeDependentSystem
import ModelingToolkit.SymbolicUtils: AbstractRule
import ModelingToolkit.Symbolics: scalarize

export RandomParameter
export RateFunction
export MarkovSystem, rate_vals, node_vals, edge_vals, rate_function,
    get_graph, steady_states, transition_matrix, generate_steady_states
export add_edge!, rem_edge!, add_vertex!, rem_vertex!, connect!
export @named, @rule

const MTK = ModelingToolkit

include("random_parameter.jl")
include("rate_functions.jl")
include("utils.jl")
include("markov_system.jl")
include("graphedits.jl")

end
