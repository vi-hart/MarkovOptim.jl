abstract type AbstractMarkovModel end

mutable struct MarkovModel{T<:Integer, P<:Real} <: AbstractMarkovModel
    graph::SimpleGraph{T}
    terms::T  # number of output terms (i.e. dimension of rate_function output)
    observed::Function
    rate_function::AbstractRateFunction  # function determining transition rates
    rate_params::RandomParameter{P}
    node_params::RandomParameter{P}
    edge_params::RandomParameter{P}
    reduced_params::Union{Nothing, Matrix{P}}
    name::Symbol
end

function MarkovModel(
    graph::SimpleGraph{T},
    terms::T,
    observed::Function,
    rate_function::RateFunction;
    rate_params::Normal{P} = Normal(),
    node_params::Normal{P} = rate_params,
    edge_params::Normal{P} = node_params,
    name::Symbol = gensym("markov")) where {T<:Integer, P<:Real}

    nd, ed = nv(graph), ne(graph)
    ratepars = RandomParameter(rate_params, terms - has_shift(rate_function), length(rate_function))
    nodepars = RandomParameter(node_params, nd, terms)
    edgepars = RandomParameter(edge_params, ed, terms)
    return MarkovModel{T, P}(graph, terms, observed, rate_function, ratepars, nodepars, edgepars, nothing, name)
end

function MarkovModel(
    graph::SimpleGraph{T};
    terms::T = one(T) << 1,
    observed::Function = (x)->x[1],
    rate_function::Symbol = :Sigmoid,
    shift::Bool = true,
    qwargs...) where T<:Integer

    ratef = getfield(MarkovOptim, rate_function)(shift)
    return MarkovModel(graph, terms, observed, ratef; qwargs...)
end

function MarkovModel(nv::T, ne::T; qwargs...) where T<:Integer
    graph = SimpleGraph{T}(nv, ne)
    return MarkovModel(graph; qwargs...)
end

function MarkovModel{T}(nv::Integer, ne::Integer; qwargs...) where T<:Integer
    graph = SimpleGraph{T}(nv, ne)
    return MarkovModel(graph; qwargs...)
end

Base.nameof(M::MarkovModel) = getfield(M, :name)
has_graph(M::MarkovModel) = hasfield(M, :graph)
get_graph(M::MarkovModel) = getfield(M, :graph)
MTK.observed(M::MarkovModel) = getfield(M, :observed)
rate_parameters(M::MarkovModel) = getfield(M, :rate_params)
node_parameters(M::MarkovModel) = getfield(M, :node_params)
edge_parameters(M::MarkovModel) = getfield(M, :edge_params)
reduced_parameters(M::MarkovModel) = getfield(M, :reduced_params)
is_reduced(M::T) where T<:MarkovModel = hasfield(T, :reduced_params) && !isnothing(M.reduced_params)

rate_function(M::MarkovModel) = getfield(M, :rate_function)
function rate_function(M::MarkovModel, x::Real)
    ps = Vector(rate_parameters(M))
    return rate_function(M)(x, ps)
end

function _reduce_params(graph::SimpleGraph, α::AbstractMatrix{P}, β::AbstractMatrix{P}) where P
    ic = _incidence(graph)
    E = ne(graph)
    D = Matrix{P}(undef, E, 2E)
    D[:,1:2:end] .= I(E)
    D[:,2:2:end] .= -I(E)
    return inv([D;abs.(D)])*[ic'α;β]
end

function reduce_params(M::MarkovModel)
    graph = get_graph(M)
    α = node_parameters(M)
    β = edge_parameters(M)
    return _reduce_params(graph, α, β)
end

function reduce_params!(M::MarkovModel)
    theta = reduce_params(M)
    setfield!(M, :reduced_params, theta)
    return M
end

function rates(M::MarkovModel, x::Real)
    !is_reduced(M) && reduce_params!(M)
    theta = reduced_parameters(M)
    return exp.(theta * rate_function(M, x))
end

function _transition_matrix(graph::SimpleGraph, rs::Vector)
    dic = _directed_incidence(graph)
    return dic*Diagonal(rs[:])*(-min.(dic, 0)')
end

transition_matrix(M::MarkovModel, x::Real) = _transition_matrix(get_graph(M), rates(M, x))

Base.Matrix(M::MarkovModel, x::Real) = transition_matrix(M, x)

function _steady_states(α::AbstractMatrix, rs::Vector)
    s = [exp(dot(rs, α[i,:])) for i in 1:size(α, 1)]
    return s /= sum(s)
end

function steady_states(M::MarkovModel, x::Real)
    rates = rate_function(M, x)
    α = node_parameters(M)
    return _steady_states(α, rates)
end

function Base.convert(
    ::Type{ODESystem},
    M::MarkovModel;
    eqs::Vector{Equation} = Equation[],
    observed::Vector{Equation} = Equation[],
    xval::Union{Nothing, Real, Function} = nothing,
    name::Symbol = nameof(M),
    kwargs...
)
#TODO: Fix equivalency with _reduce_convert
    ratepars = rate_parameters(M)
    nodepars = node_parameters(M)
    edgepars = edge_parameters(M)
    defaults = Dict()

    a,b = size(ratepars)
    @parameters γ[1:a,1:b]
    merge!(defaults, Dict(Symbolics.scalarize(γ) .=> ratepars))
    ps = [[γ[i,j] for j in Base.oneto(b)] for i in Base.oneto(a)]

    a,b = size(nodepars)
    @parameters α[1:a,1:b]
    merge!(defaults, Dict(Symbolics.scalarize(α) .=> nodepars))

    a,b = size(edgepars)
    @parameters β[1:a,1:b]
    merge!(defaults, Dict(Symbolics.scalarize(β) .=> edgepars))

    graph = get_graph(M)
    a = nv(graph)
    @variables t, x(t), s[1:a](t)
    if xval isa Real
        push!(eqs, x ~ xval)
        push!(defaults, x=>xval)
    elseif xval isa Function
        push!(eqs, x ~ xval(t))
        push!(defaults, x=>xval(0.))
    end
    D = Differential(t)

    theta = _reduce_params(graph, Symbolics.scalarize(α), Symbolics.scalarize(β))
    ratef = rate_function(M)(x, ps)
    rs = exp.(theta*ratef)

    ss = _steady_states(Symbolics.scalarize(α), ratef)
    merge!(defaults, Dict(Symbolics.scalarize(s) .=> ss))

    transition_eqs = D.(Symbolics.scalarize(s)) .~ _transition_matrix(graph, rs)*Symbolics.scalarize(s)
    append!(eqs, transition_eqs)
    push!(observed, :out ~ M.observed(s))
    sys = ODESystem(eqs, t; defaults=defaults, observed=observed, name=name, kwargs...)
    return isnothing(xval) ? sys : structural_simplify(sys)
end

function _reduce_convert(
    ::Type{ODESystem},
    M::MarkovModel;
    eqs::Vector{Equation} = Equation[],
    observed::Vector{Equation} = Equation[],
    xval::Union{Nothing, Real, Function} = nothing,
    name::Symbol = nameof(M),
    kwargs...
)
    reduce_params!(M)
    a = nv(get_graph(M))
    @variables t, x(t), s[1:a](t)
    D = Differential(t)
    _s = Symbolics.scalarize(s)
    eqs = D.(_s) .~ transition_matrix(M, x)*_s
    ss = steady_states(M, x)
    defaults = Dict(Symbolics.scalarize(s) .=> ss)
    if xval isa Real
        push!(eqs, x ~ xval)
        push!(defaults, x=>xval)
    elseif xval isa Function
        push!(eqs, x ~ xval(t))
        push!(defaults, x=>xval(0.))
    end
    push!(observed, :out ~ M.observed(s))
    sys = ODESystem(eqs, t; defaults=defaults, observed=observed, name=name, kwargs...)
    return isnothing(xval) ? sys : structural_simplify(sys)
end
