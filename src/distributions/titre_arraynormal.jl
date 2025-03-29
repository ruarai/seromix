
using Distributions
using StatsFuns

struct TitreArrayNormal{T <: Real} <: DiscreteMultivariateDistribution
    μ::AbstractVector{T}
    σ::T

    minimum::T
    maximum::T
end

function clamp_titer(x::Real, min::Real, max::Real)
    return clamp(floor(x), min, max)
end

function Distributions.rand(rng::AbstractRNG, d::TitreArrayNormal)
    return [clamp_titer(rand(rng, Normal(d.μ[i], d.σ)), d.minimum, d.maximum) for i in eachindex(d.μ)]
end

Distributions.minimum(d::TitreArrayNormal) = d.minimum
Distributions.maximum(d::TitreArrayNormal) = d.maximum
Base.length(d::TitreArrayNormal)::Int = length(d.μ)

function titre_logpdf_component(x::S, μ::T, σ::T, min::T, max::T) where {T <: Real, S <: Real}
    if x <= min
        return normlogcdf(μ, σ, x + 1)
    elseif x > min && x < max
        return log(normcdf(μ, σ, x + 1) - normcdf(μ, σ, x))
    elseif x >= max
        return normlogccdf(μ, σ, x + 1)
    else
        return typemax(T)
    end
end

function apply_logpdf(x::AbstractVector{S}, μ::SubArray{T, }, σ::T, min::T, max::T) where {T <: Real, S <: Real}
    return mapreduce(
        x -> titre_logpdf_component(x[1], x[2], σ, min, max),
        +, 
        zip(x, μ)
    )
end

function apply_logpdf(x::AbstractVector{S}, μ::Vector{T}, σ::T, min::T, max::T) where {T <: Real, S <: Real}
    return mapreduce(
        x -> titre_logpdf_component(x[1], x[2], σ, min, max),
        +, 
        zip(x, μ)
    )
end

function Distributions.logpdf(d::TitreArrayNormal{T}, x::AbstractVector{S})::T where {T <: Real, S <: Real}
    return apply_logpdf(x, d.μ, d.σ, d.minimum, d.maximum)
end