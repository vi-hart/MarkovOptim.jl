abstract type AbstractRateFunction end

struct RateFunction <: AbstractRateFunction
    f::Function
    npars::Int
    shift::Bool
end

has_shift(RF::T) where T<:RateFunction = hasfield(T, :shift) && getfield(RF, :shift)
get_function(RF::RateFunction) = getfield(RF, :f)
@inline Base.length(RF::RateFunction) = RF.npars

function (rate::RateFunction)(x::T, ps::Vector) where T
    out = map(p->get_function(rate)(x, p), ps)
    has_shift(rate) && pushfirst!(out, one(T))
    return out
end

# Simple rate functions
Sigmoid(shift::Bool=true) = RateFunction((x, p) -> tanh((x + p[1])/p[2]), 2, shift)
Linear(shift::Bool=true) = RateFunction((x, p) -> p[1]*x, 1, shift)
