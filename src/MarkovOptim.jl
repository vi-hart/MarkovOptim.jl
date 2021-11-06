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
export MarkovModel, rate_parameters, node_parameters, edge_parameters, rate_function,
    get_graph, steady_states, transition_matrix
export add_edge!, rem_edge!, add_vertex!, rem_vertex!, connect!
export @named, @rule

const MTK = ModelingToolkit

include("utils.jl")
include("random_parameter.jl")
include("rate_functions.jl")
include("markov.jl")
include("markov_system.jl")
include("graphedits.jl")

end
