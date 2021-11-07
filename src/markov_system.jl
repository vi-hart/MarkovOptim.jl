abstract type AbstractMarkovSystem <: AbstractTimeDependentSystem end

struct MarkovSystem{P<:AbstractMatrix} <: AbstractMarkovSystem
    eqs::Vector{Equation}
    """Independent variable (usually time)."""
    iv::Num
    input::Num
    output::Union{Nothing, AbstractRule}
    graph::SimpleGraph
    terms::Integer
    node_vals::P
    edge_vals::P
    rate_vals::P
    rate_func::AbstractRateFunction
    """Equations for observed variables."""
    observed::Vector{Equation}
    systems::Vector{ODESystem}
    defaults::Dict
    name::Symbol
end

function MarkovSystem(
    graph::SimpleGraph,
    terms::Integer,
    rate_func::AbstractRateFunction;
    eqs::Vector{Equation} = Equation[],
    iv::Num = only(@parameters t),
    input::Num = only(@parameters x),
    output::Union{Nothing, AbstractRule} = nothing,
    node_dist::Normal = Normal(),
    edge_dist::Normal = node_dist,
    rate_dist::Normal = edge_dist,
    observed::Vector{Equation} = Equation[],
    systems::Vector{ODESystem} = ODESystem[],
    defaults::Dict = Dict(),
    name::Symbol = gensym("markov"))

    nd, ed = nv(graph), ne(graph)
    nterms = terms - has_shift(rate_func)
    nrates = length(rate_func)
    ratevals = RandomParameter(rate_dist, nterms, nrates)
    nodevals = RandomParameter(node_dist, nd, terms)
    edgevals = RandomParameter(edge_dist, ed, terms)

    return MarkovSystem(eqs, iv, input, output, graph, terms, nodevals, edgevals, ratevals, rate_func, observed, systems, defaults, name)
end

function MarkovSystem(
    graph::SimpleGraph{T};
    terms::T = one(T) << 1,
    rate_func::Symbol = :Sigmoid,
    shift::Bool = true,
    qwargs...) where T<:Integer

    ratef = getfield(MarkovOptim, rate_func)(shift)
    return MarkovSystem(graph, terms, ratef; qwargs...)
end

function MarkovSystem(nv::T, ne::T; qwargs...) where T<:Integer
    graph = SimpleGraph{T}(nv, ne)
    return MarkovSystem(graph; qwargs...)
end

function MarkovSystem{T}(nv::Integer, ne::Integer; qwargs...) where T<:Integer
    graph = SimpleGraph{T}(nv, ne)
    return MarkovSystem(graph; qwargs...)
end

get_graph(sys::MarkovSystem) = getfield(sys, :graph)
get_input(sys::MarkovSystem) = getfield(sys, :input)
node_vals(sys::MarkovSystem) = getfield(sys, :node_vals)
edge_vals(sys::MarkovSystem) = getfield(sys, :edge_vals)
rate_vals(sys::MarkovSystem) = getfield(sys, :rate_vals)
get_rate_func(sys::MarkovSystem) = getfield(sys, :rate_func)

function node_ps(sys::MarkovSystem)
    vals = node_vals(sys)
    indices = Base.OneTo.(size(vals))
    @parameters α[indices...] = vals
    return scalarize(α)
end

function edge_ps(sys::MarkovSystem)
    vals = edge_vals(sys)
    indices = Base.OneTo.(size(vals))
    @parameters β[indices...] = vals
    return scalarize(β)
end

function rate_ps(sys::MarkovSystem)
    vals = rate_vals(sys)
    indices = Base.OneTo.(size(vals))
    @parameters γ[indices...] = vals
    return scalarize(γ)
end

function get_us(sys::MarkovSystem)
    nd = nv(get_graph(sys))
    iv = MTK.get_iv(sys)
    @variables u[1:nd](iv)
    return scalarize(u)
end

function rate_function(sys::MarkovSystem, x::Real; substitute=true)
    ps = _vecOfvec(substitute ? rate_vals(sys) : rate_ps(sys))
    return get_rate_func(sys)(x, ps)
end
rate_function(sys::MarkovSystem; kwargs...) = rate_function(sys, get_input(sys); kwargs...)

function _reduce_params(graph::SimpleGraph, α::AbstractMatrix{P}, β::AbstractMatrix{P}) where P
    ic = _incidence(graph)
    E = ne(graph)
    D = Matrix{P}(undef, E, 2E)
    D[:,1:2:end] .= I(E)
    D[:,2:2:end] .= -I(E)
    return inv([D;abs.(D)])*[ic'α;β]
end

function reduce_params(sys::MarkovSystem; substitute=true)
    graph = get_graph(sys)
    if substitute
        α = node_vals(sys)
        β = edge_vals(sys)
    else
        α = node_ps(sys)
        β = edge_ps(sys)
    end
    return _reduce_params(graph, α, β)
end

function rates(sys::MarkovSystem, x::Real; substitute=true)
    theta = reduce_params(sys; substitute=substitute)
    return exp.(theta * rate_function(sys, x; substitute=substitute))
end

function _transition_matrix(graph::SimpleGraph, rs::Vector)
    dic = _directed_incidence(graph)
    return dic*Diagonal(rs[:])*(-min.(dic, 0)')
end

function transition_matrix(sys::MarkovSystem, x::Real; substitute=true)
    return _transition_matrix(get_graph(sys), rates(sys, x; substitute=substitute))
end
function transition_matrix(sys::MarkovSystem; substitute=true)
    return _transition_matrix(get_graph(sys), rates(sys, get_input(sys); substitute=substitute))
end

function _steady_states(α::AbstractMatrix, rfs::Vector)
    s = [exp(dot(rfs, α[i,:])) for i in 1:size(α, 1)]
    return s /= sum(s)
end

function steady_states(sys::MarkovSystem, x::Real; substitute=true)
    rfs = rate_function(sys, x; substitute=substitute)
    α = substitute ? node_vals(sys) : node_ps(sys)
    return _steady_states(α, rfs)
end
steady_states(sys::MarkovSystem; substitute=true) = steady_states(sys, get_input(sys); substitute=substitute)

#TODO: Generate code for the steady states
function generate_steady_states(sys::MarkovSystem; substitute=true)
    ss = steady_states(sys; substitute=substitute)
    build_function(ss, get_input(sys))
end

function MTK.get_eqs(sys::MarkovSystem; substitute=false)
    iv = MTK.get_iv(sys)
    D = Differential(iv)
    u = get_us(sys)
    eqs = D.(u) .~ transition_matrix(sys; substitute=substitute)*u
    return append!(eqs, getfield(sys, :eqs))
end

#TODO: Include states from input variables
function MTK.get_states(sys::MarkovSystem)
    _states = Set(get_us(sys))
    for eq in getfield(sys, :eqs)
        for var in get_variables(eq.lhs)
            MTK.isparameter(var) || push!(_states, var)
        end
        for var in get_variables(eq.rhs)
            MTK.isparameter(var) || push!(_states, var)
        end
    end
    return _states
end

#TODO: return all parameters
function MTK.get_ps(sys::MarkovSystem)
    ps = vec(node_ps(sys))
    append!(ps, vec(edge_ps(sys)))
    append!(ps, vec(rate_ps(sys)))
    return ps
end

function MTK.get_observed(sys::MarkovSystem)
    out_rule = getfield(sys, :output)
    if isnothing(out_rule)
        return Equation[]
    else
        u = get_us(sys)
        iv = MTK.get_iv(sys)
        @variables out(iv)
        return out ~ out_rule(u)
    end
end
function MTK.get_defaults(sys::MarkovSystem; substitute=false)
    defaults = Dict()
    rs = rate_function(sys; substitute=substitute)
    push!(defaults, MTK.get_iv(sys)=>0.)
    merge!(defaults, getfield(sys, :defaults))
    if ~substitute
        α = node_ps(sys)
        β = edge_ps(sys)
        γ = rate_ps(sys)
        map((x, y) -> push!(defaults, x=>y), α, node_vals(sys))
        map((x, y) -> push!(defaults, x=>y), β, edge_vals(sys))
        map((x, y) -> push!(defaults, x=>y), γ, rate_vals(sys))
    end
    map((x, y) -> push!(defaults, x=>y), get_us(sys), steady_states(sys; substitute=substitute))
    return defaults
end

MTK.observed(sys::MarkovSystem) = [MTK.get_observed(sys); getfield(sys, :observed)]

MTK.calculate_jacobian(sys::MarkovSystem) = transition_matrix(sys; substitute=false)

function Base.convert(::Type{ODESystem}, sys::MarkovSystem; simplify=true, substitute=false, name=nameof(sys))
    odesys =  ODESystem(MTK.get_eqs(sys; substitute=substitute),
        MTK.get_iv(sys);
        defaults=MTK.get_defaults(sys; substitute=substitute),
        systems=MTK.get_systems(sys),
        observed=observed(sys),
        name=name
    )
    return simplify ? structural_simplify(odesys; simplify=true) : odesys
end

function Base.reduce(::Type{ODESystem}, sys::MarkovSystem; simplify=true, name=nameof(sys))
    return convert(ODESystem, sys; simplify=simplify, substitute=true, name=name)
end
