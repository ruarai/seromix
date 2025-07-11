

# A re-implementation of the adaptive Metropolis-Hastings sampler
# used in Kucharski (2018)

# It is a Metropolis-Hastings sampler with a log-normal proposal distribution
# The proposal distribution is adapted based on the acceptance rate

# In the code of Kucharski (2018), parameters are constrained during the MH sampling
# step. Here, we constrain parameters using the prior defined in the model itself.


mutable struct MHParameterSamplerOriginal <: AbstractMCMC.AbstractSampler

end

# Initial step for the sampler
# If initial parameters are provided, use them. Otherwise, 
# sample initial parameters from a uniform distribution.
function AbstractMCMC.step(
    rng::Random.AbstractRNG,
    model::AbstractMCMC.LogDensityModel,
    sampler::MHParameterSamplerOriginal;
    initial_params,
    kwargs...
)
    if !isnothing(initial_params)
        # TODO fix with link function
        theta_init = initial_params
    else
        theta_init = rand(rng, Uniform(0, 1.0), 7)
    end

    sigma_covar_0 = 0.001

    logprob_init = get_logp(theta_init, model)

    transition = ParameterSamplerTransition(theta_init, logprob_init)
    return transition, ParameterSamplerState(transition, theta_init, sigma_covar_0, 0, 0)
end

# Step function for the sampler
# Proposes a new set of parameters using a log-normal proposal distribution.
function AbstractMCMC.step(
    rng::Random.AbstractRNG,
    model::AbstractMCMC.LogDensityModel,
    sampler::MHParameterSamplerOriginal,
    state::ParameterSamplerState;
    kwargs...
)
    prev_transition = state.transition
    theta_current = prev_transition.θ

    n_accepted = state.n_accepted
    n_rejected = state.n_rejected
    n_steps = (n_accepted + n_rejected)

    target_accept_rate = 0.234
    obs_accept_rate = target_accept_rate

    # Twice n_steps to match original code
    m = n_steps * 2

    if m >= 100
        obs_accept_rate = n_accepted / n_steps
    end

    sigma_covar_adapted = max(0.00001, min(1, exp(log(state.sigma_covar) + (obs_accept_rate - target_accept_rate) * 0.999 ^ m)))

    # Get the logdensity function
    logprob_previous = get_logp(theta_current, model)

    # Theta proposal
    log_theta_current = log.(theta_current)

    log_theta_new_dist = MvNormal(log_theta_current, sigma_covar_adapted * I)
    log_theta_new = rand(rng, log_theta_new_dist)
    theta_new = exp.(log_theta_new)

    logprob_proposal = get_logp(theta_new, model)

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
function make_mh_parameter_sampler_original()
    return externalsampler(MHParameterSamplerOriginal(), adtype=Turing.DEFAULT_ADTYPE, unconstrained=false)
end

isgibbscomponent(::MHParameterSamplerOriginal) = true