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
rate_function(sys::MarkovSystem) = getfield(sys, :rate_func)
node_vals(sys::MarkovSystem) = getfield(sys, :node_vals)
edge_vals(sys::MarkovSystem) = getfield(sys, :edge_vals)
rate_vals(sys::MarkovSystem) = getfield(sys, :rate_vals)

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

function MTK.get_eqs(sys::MarkovSystem)
    graph = get_graph(sys)
    iv = MTK.get_iv(sys)
    D = Differential(iv)
    u = get_us(sys)
    x = get_input(sys)
    α = node_ps(sys)
    β = edge_ps(sys)
    γ = rate_ps(sys)
    ps = [γ[1,:] for i in 1:size(γ,1)]

    theta = _reduce_params(graph, α, β)
    rs = exp.(theta * rate_function(sys)(x, ps))
    eqs = D.(u) .~ _transition_matrix(graph, rs)*u
    return append!(eqs, getfield(sys, :eqs))
end

function reduced_eqs(sys::MarkovSystem)
    graph = get_graph(sys)
    iv = MTK.get_iv(sys)
    D = Differential(iv)
    u = get_us(sys)
    x = get_input(sys)
    α = node_vals(sys)
    β = edge_vals(sys)
    γ = rate_vals(sys)
    ps = [γ[1,:] for i in 1:size(γ,1)]

    theta = _reduce_params(graph, α, β)
    rs = exp.(theta * rate_function(sys)(x, ps))
    eqs = D.(u) .~ _transition_matrix(graph, rs)*u
    return append!(eqs, getfield(sys, :eqs))
end

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

MTK.get_ps(sys::MarkovSystem) = [node_ps(sys); edge_ps(sys); rate_ps(sys)] |> vec

function MTK.get_observed(sys::MarkovSystem)
    out_rule = getfield(sys, :output)
    if isnothing(out_rule)
        return nothing
    else
        u = get_us(sys)
        iv = MTK.get_iv(sys)
        @variables out(iv)
        return out ~ out_rule(u)
    end
end
function MTK.get_defaults(sys::MarkovSystem)
    defaults = Dict()
    x = get_input(sys)
    α = node_ps(sys)
    β = edge_ps(sys)
    γ = rate_ps(sys)
    ps = [γ[1,:] for i in 1:size(γ,1)]
    rs = rate_function(sys)(x, ps)
    push!(defaults, MTK.get_iv(sys)=>0.)
    merge!(defaults, getfield(sys, :defaults))
    map((x, y) -> push!(defaults, x=>y), α, node_vals(sys))
    map((x, y) -> push!(defaults, x=>y), β, edge_vals(sys))
    map((x, y) -> push!(defaults, x=>y), γ, rate_vals(sys))
    map((x, y) -> push!(defaults, x=>y), get_us(sys), _steady_states(α, rs))
    return defaults
end

function reduced_defaults(sys::MarkovSystem)
    defaults = Dict()
    x = get_input(sys)
    α = node_vals(sys)
    β = edge_vals(sys)
    γ = rate_vals(sys)
    ps = [γ[1,:] for i in 1:size(γ,1)]
    rs = rate_function(sys)(x, ps)
    push!(defaults, MTK.get_iv(sys)=>0.)
    merge!(defaults, getfield(sys, :defaults))
    map((x, y) -> push!(defaults, x=>y), get_us(sys), _steady_states(α, rs))
    return defaults
end
# function MTK.equations(sys::MarkovSystem) end
# function MTK.states(sys::MarkovSystem) end
# function MTK.parameters(sys::MarkovSystem) end
MTK.observed(sys::MarkovSystem) = [MTK.get_observed(sys); getfield(sys, :observed)]

function Base.convert(::Type{ODESystem}, sys::MarkovSystem; simplify=true, name=nameof(sys))
    odesys =  ODESystem(MTK.get_eqs(sys),
        MTK.get_iv(sys);
        defaults=MTK.get_defaults(sys),
        systems=MTK.get_systems(sys),
        observed=observed(sys),
        name=name
    )
    return simplify ? structural_simplify(odesys; simplify=true) : odesys
end

function Base.reduce(::Type{ODESystem}, sys::MarkovSystem; simplify=true, name=nameof(sys))
    odesys = ODESystem(reduced_eqs(sys),
        MTK.get_iv(sys);
        defaults=reduced_defaults(sys),
        systems=MTK.get_systems(sys),
        observed=observed(sys),
        name=name
    )
    return simplify ? structural_simplify(odesys; simplify=true) : odesys
end
