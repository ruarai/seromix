mutable struct MHInfectionSampler <: AbstractMCMC.AbstractSampler
end

import Turing.Inference: isgibbscomponent

isgibbscomponent(::MHInfectionSampler) = true

struct Transition{T}
    θ::AbstractVector{T}
end

struct SamplerState{P<:Transition, T}
    transition::P
    θ::AbstractVector{T}
end

AbstractMCMC.getparams(state::SamplerState) = state.θ
function AbstractMCMC.setparams!!(state::SamplerState, θ)
    return SamplerState(
        Transition(θ), θ
    )
end


function AbstractMCMC.step(
    rng::Random.AbstractRNG,
    model::AbstractMCMC.LogDensityModel,
    spl::MHInfectionSampler;
    kwargs...
)
    d = LogDensityProblems.dimension(model.logdensity)

    theta_init = convert(Vector{Real}, fill(false, d))

    transition = Transition(theta_init)
    return transition, SamplerState(transition, theta_init)
end

function AbstractMCMC.step(
    rng::Random.AbstractRNG,
    model::AbstractMCMC.LogDensityModel,
    sampler::MHInfectionSampler,
    state::SamplerState;
    kwargs...
)
    prev_transition = state.transition
    theta = prev_transition.θ

    theta_proposal = propose_theta(rng, theta, model)

    logprob_previous = LogDensityProblems.logdensity(model.logdensity, theta)
    logprob_proposal = LogDensityProblems.logdensity(model.logdensity, theta_proposal)

    # logratio_proposal = 0.0

    log_α = logprob_proposal - logprob_previous# + logratio_proposal

    if -Random.randexp(rng) <= log_α
        transition = Transition(theta_proposal)
        return transition, SamplerState(transition, theta_proposal)
    else
        return prev_transition, SamplerState(prev_transition, theta)
    end
end

# function propose_theta(rng, theta, model)
#     p_swap = rand(rng, Beta(2, length(theta) - 2))

#     theta_mask = rand(rng, Bernoulli(p_swap), length(theta))
#     theta_prop = xor.(theta, theta_mask)

#     return Vector{Real}(theta_prop)
# end

function propose_theta(rng, theta::AbstractVector{T}, model) where T <: Real
    n = length(theta)
    p_swap = rand(rng, Beta(2, n - 2))
    
    # Generate mask directly (single allocation)
    theta_mask = rand(rng, Bernoulli(p_swap), n)
    
    # Use pre-allocated result vector to avoid intermediate allocation
    theta_prop = similar(theta)
    
    # Apply XOR operation directly
    @inbounds for i in 1:n
        theta_prop[i] = xor(theta[i], theta_mask[i])
    end
    
    return Vector{Real}(theta_prop)
end