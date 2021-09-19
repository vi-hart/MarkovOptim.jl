abstract type AbstractMarkovModel end

mutable struct MarkovModel{T<:Integer, P<:Real} <: AbstractMarkovModel
    graph::SimpleGraph{T}
    open::T  # number of open states
    terms::T  # number of output terms (i.e. dimension of rate_function output)
    rate_function::AbstractRateFunction  # function determining transition rates
    rate_params::RandomParameter{P}
    node_params::RandomParameter{P}
    edge_params::RandomParameter{P}
end

function MarkovModel(
    graph::SimpleGraph{T},
    open::T,
    terms::T,
    rate_function::RateFunction;
    rate_params::Normal{P} = Normal(),
    node_params::Normal{P} = rate_params,
    edge_params::Normal{P} = node_params) where {T<:Integer, P<:Real}

    nd, ed = nv(graph), ne(graph)
    ratepars = RandomParameter(rate_params, terms - has_shift(rate_function), length(rate_function))
    nodepars = RandomParameter(node_params, nd, terms)
    edgepars = RandomParameter(edge_params, ed, terms)
    return MarkovModel{T, P}(graph, open, terms, rate_function, ratepars, nodepars, edgepars)
end

function MarkovModel(
    graph::SimpleGraph{T};
    open::T = one(T),
    terms::T = one(T) << 1,
    rate_function::Symbol = :sigmoid,
    qwargs...) where T<:Integer

    ratef = getfield(MarkovOptim, rate_function)
    return MarkovModel(graph, open, terms, ratef; qwargs...)
end

function MarkovModel(nv::T, ne::T; qwargs...) where T<:Integer
    graph = SimpleGraph{T}(nv, ne)
    return MarkovModel(graph; qwargs...)
end

function MarkovModel{T}(nv::Integer, ne::Integer; qwargs...) where T<:Integer
    graph = SimpleGraph{T}(nv, ne)
    return MarkovModel(graph; qwargs...)
end

has_graph(M::MarkovModel) = hasfield(M, :graph)
get_graph(M::MarkovModel) = getfield(M, :graph)
rate_parameters(M::MarkovModel) = getfield(M, :rate_params)
node_parameters(M::MarkovModel) = getfield(M, :node_params)
edge_parameters(M::MarkovModel) = getfield(M, :edge_params)
