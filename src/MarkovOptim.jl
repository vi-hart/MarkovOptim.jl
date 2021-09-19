module MarkovOptim

using Distributions

export RandomParameter
export RateFunction

include("random_parameter.jl")
include("rate_functions.jl")

end
