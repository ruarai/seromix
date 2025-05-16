mutable struct MHParameterSampler <: AbstractMCMC.AbstractSampler

end


isgibbscomponent(::MHParameterSampler) = true

struct ParameterSamplerTransition{T}
    θ::AbstractVector{T}
end

struct ParameterSamplerState{P<:ParameterSamplerTransition, T}
    transition::P
    θ::AbstractVector{T}

    sigma_covar::Float64
    n_accepted::Int
    n_rejected::Int
end

AbstractMCMC.getparams(state::ParameterSamplerState) = state.θ
function AbstractMCMC.setparams!!(state::ParameterSamplerState, θ)
    return ParameterSamplerState(ParameterSamplerTransition(θ), θ, state.sigma_covar, state.n_accepted, state.n_rejected)
end


function AbstractMCMC.step(
    rng::Random.AbstractRNG,
    model::AbstractMCMC.LogDensityModel,
    sampler::MHParameterSampler;
    kwargs...
)
    d = LogDensityProblems.dimension(model.logdensity)

    theta_init = [2.0, 2.0, 0.15, 0.05, 0.05] .* rand(rng, Uniform(0.95, 1.05), 5)

    sigma_covar_0 = 0.001

    println(theta_init)

    transition = ParameterSamplerTransition(theta_init)
    return transition, ParameterSamplerState(transition, theta_init, sigma_covar_0, 0, 0)
end


function AbstractMCMC.step(
    rng::Random.AbstractRNG,
    model::AbstractMCMC.LogDensityModel,
    sampler::MHParameterSampler,
    state::ParameterSamplerState;
    kwargs...
)
    prev_transition = state.transition
    theta_current = prev_transition.θ

    n_accepted = state.n_accepted
    n_rejected = state.n_rejected
    n_steps = n_accepted + n_rejected

    target_accept_rate = 0.234
    obs_accept_rate = target_accept_rate

    if n_steps >= 100
        obs_accept_rate = n_accepted / n_steps
    end

    sigma_covar_adapted = max(0.00001, min(1, exp(log(state.sigma_covar) + (obs_accept_rate - target_accept_rate) * 0.999 ^ n_steps)))

    # Get the logdensity function
    f = model.logdensity

    context = DynamicPPL.DefaultContext()

    varinfo_prev = DynamicPPL.unflatten(f.varinfo, theta_current)
    logprob_previous = DynamicPPL.getlogp(last(DynamicPPL.evaluate!!(f.model, varinfo_prev, context)))

    # Theta proposal
    log_theta_current = log.(theta_current)

    log_theta_new_dist = MvNormal(log_theta_current, sigma_covar_adapted)
    log_theta_new = rand(rng, log_theta_new_dist)
    theta_new = exp.(log_theta_new)


    # TODO replace with direct function call?
    varinfo_proposal = DynamicPPL.unflatten(f.varinfo, theta_new)
    logprob_proposal = DynamicPPL.getlogp(last(DynamicPPL.evaluate!!(f.model, varinfo_proposal, context)))


    log_target_ratio = logprob_proposal - logprob_previous

    # Calculate log of Hastings correction term: log(H)
    # H = q(theta_current | theta_new) / q(theta_new | theta_current)
    # For log-normal proposal: H = (PRODUCT theta_new_i) / (PRODUCT theta_current_i)
    # log(H) = sum(log.(theta_new)) - sum(log.(theta_current))
    # which is sum(log_theta_new) - sum(log_theta_current)

    # NOTE: this may only be valid where priors are all uniform?

    log_hastings_ratio = sum(log_theta_new) - sum(log_theta_current)
    acceptance_ratio = log_target_ratio + log_hastings_ratio

    if -Random.randexp(rng) <= acceptance_ratio
        # Accept theta_new
        transition = ParameterSamplerTransition(theta_new)
        n_accepted += 1
        return transition, ParameterSamplerState(transition, theta_new, sigma_covar_adapted, n_accepted, n_rejected)
    else
        # Reject theta_new, return theta_current
        transition = ParameterSamplerTransition(theta_current)
        n_rejected += 1
        return transition, ParameterSamplerState(transition, theta_current, sigma_covar_adapted, n_accepted, n_rejected)
    end
end

function make_mh_parameter_sampler()
    return externalsampler(MHParameterSampler(), adtype=Turing.DEFAULT_ADTYPE, unconstrained=false)
end