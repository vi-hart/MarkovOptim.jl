module MarkovOptim

using Distributions
using LightGraphs
using LinearAlgebra

export RandomParameter
export RateFunction
export MarkovModel

include("utils.jl")
include("random_parameter.jl")
include("rate_functions.jl")
include("markov.jl")

end
