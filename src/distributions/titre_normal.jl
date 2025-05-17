
using Distributions
using StatsFuns

struct TitreNormal{T <: Real} <: Distribution{Univariate, Discrete}
    μ::T
    σ::T

    minimum::T
    maximum::T
end

function Distributions.rand(rng::AbstractRNG, d::TitreNormal)
    return clamp_titer(rand(rng, Normal(d.μ, d.σ)), d.minimum, d.maximum)
end

Distributions.minimum(d::TitreNormal) = d.minimum
Distributions.maximum(d::TitreNormal) = d.maximum
Base.length(d::TitreNormal)::Int = length(d.μ)

function titre_logpdf_component(x::S, μ::T, σ::T, min::T, max::T) where {T <: Real, S <: Real}
    if x <= min
        return normlogcdf(μ, σ, x + 1)
    elseif x < max
        return logsubexp(normlogcdf(μ, σ, x + 1), normlogcdf(μ, σ, x))
    elseif x >= max
        return normlogccdf(μ, σ, x)
    else
        return typemax(T)
    end
end

function Distributions.logpdf(d::TitreNormal{T}, x::S)::T where {T <: Real, S <: Real}
    return titre_logpdf_component(x, d.μ, d.σ, d.minimum, d.maximum)
end