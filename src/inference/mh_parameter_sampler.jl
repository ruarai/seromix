

# A re-implementation of the adaptive Metropolis-Hastings sampler
# used in Kucharski (2018)

# It is a Metropolis-Hastings sampler with a log-normal proposal distribution
# The proposal distribution is adapted based on the acceptance rate

# In the code of Kucharski (2018), parameters are constrained during the MH sampling
# step. Here, we constrain parameters using the prior defined in the model itself.


mutable struct MHParameterSampler <: AbstractMCMC.AbstractSampler

end

struct ParameterSamplerTransition{T, L <: Real}
    θ::AbstractVector{T}
    lp::L
end

struct ParameterSamplerState{P<:ParameterSamplerTransition, T}
    transition::P
    θ::AbstractVector{T}

    sigma_covar::Float64
    n_accepted::Int
    n_rejected::Int
end


# Initial step for the sampler
# If initial parameters are provided, use them. Otherwise, 
# sample initial parameters from a uniform distribution.
function AbstractMCMC.step(
    rng::Random.AbstractRNG,
    model::AbstractMCMC.LogDensityModel,
    sampler::MHParameterSampler;
    initial_params,
    kwargs...
)
    if !isnothing(initial_params)
        theta_init = initial_params
    else
        theta_init = rand(rng, Uniform(0, 1.0), 7)
    end

    sigma_covar_0 = 0.001

    varinfo_init = DynamicPPL.unflatten(model.logdensity.varinfo, theta_init)
    logprob_init = DynamicPPL.getlogp(last(DynamicPPL.evaluate!!(model.logdensity.model, varinfo_init, DynamicPPL.DefaultContext())))

    transition = ParameterSamplerTransition(theta_init, logprob_init)
    return transition, ParameterSamplerState(transition, theta_init, sigma_covar_0, 0, 0)
end

# Step function for the sampler
# Proposes a new set of parameters using a log-normal proposal distribution.
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

    # Correction term for log-normal proposals around theta
    log_hastings_ratio = sum(log_theta_new) - sum(log_theta_current)
    acceptance_ratio = log_target_ratio + log_hastings_ratio

    if -Random.randexp(rng) <= acceptance_ratio
        # Accept theta_new
        transition = ParameterSamplerTransition(theta_new, logprob_proposal)
        n_accepted += 1
        return transition, ParameterSamplerState(transition, theta_new, sigma_covar_adapted, n_accepted, n_rejected)
    else
        # Reject theta_new, return theta_current
        transition = ParameterSamplerTransition(theta_current, logprob_previous)
        n_rejected += 1
        return transition, ParameterSamplerState(transition, theta_current, sigma_covar_adapted, n_accepted, n_rejected)
    end
end


# Below is necessary framework for use in Turing

# Puts the sampler into an external sampler for use in Turing
function make_mh_parameter_sampler()
    return externalsampler(MHParameterSampler(), adtype=Turing.DEFAULT_ADTYPE, unconstrained=false)
end

isgibbscomponent(::MHParameterSampler) = true

AbstractMCMC.getparams(state::ParameterSamplerState) = state.θ
function AbstractMCMC.setparams!!(state::ParameterSamplerState, θ)
    return ParameterSamplerState(ParameterSamplerTransition(θ, state.transition.lp), θ, state.sigma_covar, state.n_accepted, state.n_rejected)
end