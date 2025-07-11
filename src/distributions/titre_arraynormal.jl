

struct TitreArrayNormal{T <: Real} <: DiscreteMultivariateDistribution
    μ::AbstractVector{T}
    σ::T

    minimum::T
    maximum::T

    corrected::Bool
end

TitreArrayNormal(μ, σ, minimum, maximum) = TitreArrayNormal(μ, σ, minimum, maximum, true)

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
        return normlogcdf(μ, σ, min + 1) # P(Y <= 1) for x <= min
    elseif x >= max
        return normlogccdf(μ, σ, max) # P(Y > 8) for x >= 8
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

function apply_logpdf_simd(x, μ, σ, min_val, max_val)
    l_sum = 0.0
    inv_sigma_sqrt2 = 1.0 / (σ * sqrt(2.0))
    @turbo for i in eachindex(x)
        arg1 = (x[i] - μ[i]) * inv_sigma_sqrt2
        arg2 = (x[i] + 1.0 - μ[i]) * inv_sigma_sqrt2
        arg_min = (min_val + 1.0 - μ[i]) * inv_sigma_sqrt2
        arg_max = (max_val - μ[i]) * inv_sigma_sqrt2

        # Branch 1: x <= min
        # result1 = normlogcdf(μ, σ, min + 1.0)
        logcdf_min = log(0.5 * (1.0 + erf(arg_min)))

        # Branch 2: x >= max
        # result2 = normlogccdf(μ, σ, max)
        logccdf_max = log(0.5 * (1.0 - erf(arg_max)))

        # Branch 3: min < x < max
        # result3 = logsubexp(normlogcdf(μ, σ, x + 1.0), normlogcdf(μ, σ, x))
        logcdf_x1 = log(0.5 * (1.0 + erf(arg2)))
        logcdf_x = log(0.5 * (1.0 + erf(arg1)))
        logsubexp_val = logcdf_x1 + log1p(-exp(logcdf_x - logcdf_x1))

        mask1 = x[i] <= min_val
        mask2 = x[i] >= max_val

        l_sum += ifelse(mask1, logcdf_min, ifelse(mask2, logccdf_max, logsubexp_val))
    end
    return ifelse(isnan(l_sum), -Inf, l_sum)
end

# A more direct translation of Kucharski
# note the treatment of max
function titre_logpdf_component_uncorrected(x::S, μ::T, σ::T, min, max) where {T <: Real, S <: Real}
    if x <= min
        return normlogcdf(μ, σ, min + 1) # P(Y <= 1) for x <= min
    elseif x > max
        return normlogccdf(μ, σ, max) # P(Y > 8) for x > 8 (i.e. x => 9)
    else
        return logsubexp(normlogcdf(μ, σ, x + 1), normlogcdf(μ, σ, x)) # P(x < Y <= x+1)
    end
end

function apply_logpdf_uncorrected(x, μ, σ, min, max)
    l_sum = 0
    @inbounds for i in eachindex(x)
        l_sum += titre_logpdf_component_uncorrected(x[i], μ[i], σ, min, max)
    end
    return l_sum
end


function Distributions.logpdf(d::TitreArrayNormal{T}, x::AbstractVector{S})::T where {T <: Real, S <: Real}
    if d.corrected
        return apply_logpdf_simd(x, d.μ, d.σ, d.minimum, d.maximum)
    else
        return apply_logpdf_uncorrected(x, d.μ, d.σ, d.minimum, d.maximum)
    end
end