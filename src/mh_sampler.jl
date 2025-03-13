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

    logratio_proposal = 0.0

    log_α = logprob_proposal - logprob_previous + logratio_proposal

    if -Random.randexp(rng) <= log_α
        transition = Transition(theta_proposal)
        return transition, SamplerState(transition, theta_proposal)
    else
        return prev_transition, SamplerState(prev_transition, theta)
    end
end

function propose_theta(rng, theta, model)

    # How to get this info here?
    n_strain = 3
    n_t_steps = 50

    theta_mat = reshape(theta, n_t_steps, n_strain)
    theta_mat_prop = copy(theta_mat)

    ix_strain = sample(rng, 1:n_strain)

    if any(theta_mat[:, ix_strain])
        theta_mat_prop[:, ix_strain] .= false

        if rand() < 0.5
            ix_t = sample(rng, 1:n_t_steps)
            theta_mat_prop[ix_t, ix_strain] = true
        end

    else
        ix_t = sample(rng, 1:n_t_steps)
        theta_mat_prop[ix_t, ix_strain] = true
    end

    return Vector{Real}(vec(theta_mat_prop))
end