module MarkovOptim

using Distributions
using LightGraphs

export RandomParameter
export RateFunction
export MarkovModel

include("random_parameter.jl")
include("rate_functions.jl")
include("markov.jl")

end
