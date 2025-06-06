

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


# A more direct translation of Kucharski
# Note that this treats max a bit weird.
function titre_logpdf_component(x::S, μ::T, σ::T, min::T, max::T) where {T <: Real, S <: Real}
    if x <= min
        return normlogcdf(μ, σ, min + 1) # P(Y <= 1) for x <= min
    elseif x > max
        return normlogccdf(μ, σ, max) # P(Y > 8) for x > 8
    else
        return logsubexp(normlogcdf(μ, σ, x + 1), normlogcdf(μ, σ, x)) # P(x < Y <= x+1)
    end
end


function apply_logpdf(x, μ, σ, min, max)
    l_sum = 0

    @inbounds for i in eachindex(x)
        l_sum += titre_logpdf_component(x[i], μ[i], σ, min, max)
    end

    return l_sum
end

function Distributions.logpdf(d::TitreArrayNormal{T}, x::AbstractVector{S})::T where {T <: Real, S <: Real}
    return apply_logpdf(x, d.μ, d.σ, d.minimum, d.maximum)
end