

mutable struct MHInfectionSampler <: AbstractMCMC.AbstractSampler

end

struct Transition{T}
    θ::AbstractVector{T}
end

struct SamplerState{T<:Transition}
    "Current [`Transition`](@ref)."
    transition::T
end

function AbstractMCMC.step(
    rng::Random.AbstractRNG,
    model::AbstractMCMC.LogDensityModel,
    spl::MHInfectionSampler;
    kwargs...
)
    # init_theta = [randn(rng), false]
    theta_init = [0]

    d = LogDensityProblems.dimension(model.logdensity)


    println("A, d = $d")
    println("theta_init = $theta_init")

    transition = Transition(theta_init)
    state = SamplerState(transition)


    println("transition = $transition")

    return transition, state
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

    log_α = logprob_proposal - logprob_previous

    println("B")

    if -Random.randexp(rng) <= log_α
        transition = Transition(theta_proposal)
        return transition, SamplerState(transition)
    else
        return prev_transition, SamplerState(prev_transition)
    end
end


function propose_theta(rng, theta, model)
    d = LogDensityProblems.dimension(model.logdensity)

    theta_proposal = theta

    println("proposal = $theta_proposal")

    # theta_proposal[2] = rand(rng, Bernoulli(0.5))

    return theta_proposal
end