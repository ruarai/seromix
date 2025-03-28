
using Distributions

struct TitreNormal{T, S} <: DiscreteUnivariateDistribution
    normal_dist::T

    minimum::S
    maximum::S
end

function Distributions.rand(rng::AbstractRNG, d::TitreNormal)
    y = rand(rng, d.normal_dist)

    return Int(clamp(floor(y), d.minimum, d.maximum))
end

Distributions.minimum(d::TitreNormal) = d.minimum
Distributions.maximum(d::TitreNormal) = d.maximum

function Distributions.insupport(d::TitreNormal, x::Real)
    return d.minimum <= x <= d.maximum
end



function Distributions.logpdf(d::TitreNormal, x::Real)
    return(log(pdf(d, x))) # Could be optimised?
end

function Distributions.pdf(d::TitreNormal, x::Real)
    if x == d.minimum
        return(cdf(d.normal_dist, x + 1))
    elseif x > d.minimum && x < d.maximum
        return(cdf(d.normal_dist, x + 1) - cdf(d.normal_dist, x))
    elseif x == d.maximum
        return(1 - cdf(d.normal_dist, x))
    end

    return zero(x)
end