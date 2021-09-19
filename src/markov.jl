abstract type AbstractMarkovModel end

mutable struct MarkovModel{T<:Integer, P<:Real} <: AbstractMarkovModel
    graph::SimpleGraph{T}
    open::T  # number of open states
    terms::T  # number of output terms (i.e. dimension of rate_function output)
    rate_function::AbstractRateFunction  # function determining transition rates
    rate_params::RandomParameter{P}
    node_params::RandomParameter{P}
    edge_params::RandomParameter{P}
    reduced_params::Union{Nothing, Matrix{P}}
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
    return MarkovModel{T, P}(graph, open, terms, rate_function, ratepars, nodepars, edgepars, nothing)
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
reduced_parameters(M::MarkovModel) = getfield(M, :reduced_params)
is_reduced(M::T) where T<:MarkovModel = hasfield(T, :reduced_params) && !isnothing(M.reduced_params)

rate_function(M::MarkovModel) = getfield(M, :rate_function)
function rate_function(M::MarkovModel{T, P}, x::P) where {T, P}
    ps = get_parameters(rate_parameters(M))
    return rate_function(M)(x, ps)
end

function reduce_params(M::MarkovModel{T, P}) where {T, P}
    graph = get_graph(M)
    ic = _incidence(graph)
    α = node_parameters(M)
    β = edge_parameters(M)
    E = ne(graph)
    D = Matrix{P}(undef, E, 2E)
    D[:,1:2:end] .= I(E)
    D[:,2:2:end] .= -I(E)
    return inv([D;abs.(D)])*[ic'α;β]
end

function reduce_params!(M::MarkovModel)
    theta = reduce_params(M)
    setfield!(M, :reduced_params, theta)
    return M
end

function rates(M::MarkovModel{T, P}, x::P) where {T, P}
    !is_reduced(M) && reduce_params!(M)
    theta = reduced_parameters(M)
    return exp.(theta * rate_function(M, x))
end

function transition_matrix(M::MarkovModel{T, P}, x::P) where {T, P}
    r = rates(M, x)
    dic = _directed_incidence(get_graph(M))
    return dic*Diagonal(r[:])*(-min.(dic, 0)')
end
Base.Matrix(M::MarkovModel{T, P}, x::P) where {T, P} = transition_matrix(M, x)

function steady_states(M::MarkovModel{T, P}, x::P) where {T, P}
    rates = rate_function(M, x)
    α = node_parameters(M)
    s = [exp(dot(rates, α[i,:])) for i in 1:size(α, 1)]
    return s /= sum(s)
end
