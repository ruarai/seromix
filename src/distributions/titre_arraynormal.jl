
using Distributions
using SpecialFunctions

struct TitreArrayNormal{T <: Real} <: DiscreteMultivariateDistribution
    μ::AbstractVector{T} # May cause issues with type stability?
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

# function Distributions.insupport(d::TitreMvNormal, x::Real)
#     return d.minimum <= x <= d.maximum
# end


function titre_logpdf_component(x::S, μ::T, σ::T, min::T, max::T) where {T <: Real, S <: Real}
    if x <= min
        return log(normal_cdf(x + 1, μ, σ))
    elseif x > min && x < max
        return log(normal_cdf(x + 1, μ, σ) - normal_cdf(x, μ, σ))
    elseif x >= max
        return log(1 - normal_cdf(x, μ, σ))
    else
        return typemax(T)
    end
end

function normal_cdf(x::S, μ::T, σ::T) where {T <: Real, S <:Real}
    return 0.5 * (1.0 + SpecialFunctions.erf((x - μ) / (σ * 1.4142135623730951)))
end

function apply_logcdf(x::AbstractArray{S, M}, μ::SubArray{T, }, σ::T, min::T, max::T) where {M, T <: Real, S <: Real}
    lp::T = 0.0

    @inbounds for i::Int in 1:length(μ)
        lp += titre_logpdf_component(x[i], μ[i], σ, min, max)
    end

    return lp
end

function apply_logcdf(x::AbstractArray{S, M}, μ::Vector{T}, σ::T, min::T, max::T) where {M, T <: Real, S <: Real}
    lp::T = 0.0

    @inbounds for i::Int in 1:length(μ)
        lp += titre_logpdf_component(x[i], μ[i], σ, min, max)
    end

    return lp
end

function Distributions.logpdf(d::TitreArrayNormal{T}, x::AbstractArray{S, M})::T where {M, T <: Real, S <: Real}
    return apply_logcdf(x, d.μ, d.σ, d.minimum, d.maximum)
end


# a = TitreArrayNormal(fill(3.0, 100), 1.5, 0.0, 8.0)

# y = rand(a)
# logpdf(a, y)

# @code_warntype logpdf(a, y)

# using BenchmarkTools

# @benchmark logpdf(a, y) setup=(a = a, y = rand(a))

# @benchmark logpdf(a, y) setup=(a = a, y = view(rand(a), :))

# rand(a)


# @model function test_model(y)

#     a ~ Normal(0, 2)

#     y ~ TitreArrayNormal(fill(a, length(y)), 1.5, const_titre_min, const_titre_max)
# end

# model = test_model(y)

# sample(model, MH(), 100)