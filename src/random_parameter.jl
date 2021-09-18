mutable struct RandomParameter{T} <: AbstractMatrix{T}
    ps::Vector{Vector{T}}
    dist::Normal{T}
end

function RandomParameter(vec::Vector{VT}, dist::Normal{T}) where {T, VT<:Vector{T}}
    return RandomParameter{T}(vec, dist)
end

RandomParameter(::Type{T}) where T = RandomParameter(Normal(zero(T), one(T)), 0, 0)
RandomParameter() = RandomParameter(Float64)
RandomParameter(::Type{T}, n::Int, m::Int) where T = RandomParameter(Normal(zero(T), one(T)), n, m)
RandomParameter(n::Int, m::Int) = RandomParameter(Float64, n, m)
function RandomParameter(dist::Normal, n::Int, m::Int)
    vec = map(x->rand(dist, m), 1:n)
    RandomParameter(vec, dist)
end

has_parameters(RP::T) where T<:RandomParameter = hasfield(T, :ps)
get_parameters(RP::RandomParameter) = getfield(RP, :ps)
has_distribution(RP::T) where T<:RandomParameter = hasfield(T, :dist)
get_distribution(RP::RandomParameter) = getfield(RP, :dist)
Distributions.Distribution(RP::RandomParameter) = getfield(RP, :dist)
Distributions.params(RP::RandomParameter) = Distributions.params(getfield(RP, :dist))

@inline Base.size(RP::RandomParameter) = (length(RP.ps), size(RP.ps[1])...)
@inline Base.length(RP::RandomParameter) = prod(size(RP))
@inline Base.eachindex(RP::RandomParameter) = Base.OneTo(length(RP))
@inline Base.IteratorSize(RP::RandomParameter) = Base.HasLength()

Base.@propagate_inbounds Base.getindex(RP::RandomParameter{T}, I1::Int, I2::Int) where {T} = RP.ps[I1][I2]
Base.@propagate_inbounds function Base.getindex(RP::RandomParameter{T}, I::Int) where {T}
    return getindex(RP, fldmod1(I, size(RP, 2))...)
end
Base.@propagate_inbounds Base.setindex!(RP::RandomParameter{T}, v, I1::Int, I2::Int) where {T} = RP.ps[I1][I2] = v
Base.@propagate_inbounds function Base.setindex!(RP::RandomParameter{T}, v, I::Int) where {T}
    return setindex(RP, fldmod1(I, size(RP, 2))...)
end

recursivecopy(RP::RandomParameter) = RandomParameter(copy.(RP.ps), deepcopy(RP.dist))

Base.rand(RP::RandomParameter{T}) where T = rand(RP.dist, size(RP, 2))
Base.rand(RP::RandomParameter{T}, I::Integer) where T = map(i->rand(RP), Base.OneTo(I))

@inline Base.push!(RP::RandomParameter{T}) where T = push!(RP, rand(RP))
function Base.push!(RP::RandomParameter{T}, item::Vector) where T
    if size(RP, 2) != length(item)
        throw(ArgumentError("number of rows of each array must match (got $((size(RP, 2), length(item))))"))
    end
    push!(RP.ps, item)
    return RP
end

@inline Base.insert!(RP::RandomParameter{T}, I::Integer) where T = insert!(RP, I, rand(RP))
function Base.insert!(RP::RandomParameter{T}, I::Integer, item::Vector) where T
    if size(RP, 2) != length(item)
        throw(ArgumentError("number of rows of each array must match (got $((size(RP, 2), length(item))))"))
    end
    insert!(RP.ps, I, item)
    return RP
end

@inline Base.append!(RP::RandomParameter{T}, I::Integer) where T = append!(RP, rand(RP, I))
function Base.append!(RP::RandomParameter{T}, items::Vector) where T
    for item in copy(items)
        push!(RP, item)
    end
    return RP
end

function Base.pop!(RP::RandomParameter{T}) where T
    pop!(RP.ps)
    return RP
end

function Base.popat!(RP::RandomParameter{T}, I::Integer) where T
    popat!(RP.ps, I)
    return RP
end
