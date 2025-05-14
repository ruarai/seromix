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

    theta_init = [4.0, 2.0, 0.15, 0.05, 0.05] .* rand(rng, Uniform(0.95, 1.05), 5)

    sigma_covar_0 = 0.001

    println(theta_init)

    transition = ParameterSamplerTransition(theta_init)
    return transition, ParameterSamplerState(transition, theta_init, sigma_covar_0, 0, 0)
end




# # Adaptive covariance matrix
# if(m==1){
# }else{
#   epsilon0=max(0.00001,min(1,exp(log(epsilon0)+(accept_rateT-0.234)*0.999^m)))
#   #cov_matrix_theta=epsilon0*cov_matrix_thetaA
#   cov_matrix_basic=epsilon0*cov_matrix_theta0
#   #varpart_prob0=max(0.02,min(0.25,exp(log(varpart_prob0)+(accept_rateH-0.234)*0.999^m))) # resample max of 25%, min of 2%
# }
#   theta_star = as.numeric(exp(mvrnorm(1,log(theta_initial), Sigma=covarbasic)))
  
#   names(theta_star)=names(theta_initial)
  
#   # Check parameters biologically plausible
#   if(theta_star[["wane"]]>1 & sum(pmask=="wane")==0 ){ likelihoodOK = 0 } # theta_star[["error"]]<0 | theta_star[["tau1"]]<0 | theta_star[["tau2"]]<0 theta_star[["mu"]]<0 | theta_star[["mu"]]>10 |theta_star[["muShort"]]<0 | theta_star[["sigma"]]<0 | theta_star[["sigma2"]]<0  | theta_star[["wane"]]<0 | 
#   if(theta_star[["mu"]]>10 & sum(pmask=="mu")==0 ){ likelihoodOK = 0 } # theta_star[["error"]]<0 | theta_star[["tau1"]]<0 | theta_star[["tau2"]]<0 theta_star[["mu"]]<0 | theta_star[["mu"]]>10 |theta_star[["muShort"]]<0 | theta_star[["sigma"]]<0 | theta_star[["sigma2"]]<0  | theta_star[["wane"]]<0 | 
#   if(theta_star[["muShort"]]>10 & sum(pmask=="muShort")==0 ){ likelihoodOK = 0 } # theta_star[["error"]]<0 | theta_star[["tau1"]]<0 | theta_star[["tau2"]]<0 theta_star[["mu"]]<0 | theta_star[["mu"]]>10 |theta_star[["muShort"]]<0 | theta_star[["sigma"]]<0 | theta_star[["sigma2"]]<0  | theta_star[["wane"]]<0 | 
#   if(theta_star[["error"]]>10 & sum(pmask=="error")==0 ){ likelihoodOK = 0 } # theta_star[["error"]]<0 | theta_star[["tau1"]]<0 | theta_star[["tau2"]]<0 theta_star[["mu"]]<0 | theta_star[["mu"]]>10 |theta_star[["muShort"]]<0 | theta_star[["sigma"]]<0 | theta_star[["sigma2"]]<0  | theta_star[["wane"]]<0 | 
#   if(theta_star[["sigma"]]>10 & sum(pmask=="sigma")==0 ){ likelihoodOK = 0 } # theta_star[["error"]]<0 | theta_star[["tau1"]]<0 | theta_star[["tau2"]]<0 theta_star[["mu"]]<0 | theta_star[["mu"]]>10 |theta_star[["muShort"]]<0 | theta_star[["sigma"]]<0 | theta_star[["sigma2"]]<0  | theta_star[["wane"]]<0 | 
#   if(theta_star[["sigma2"]]>10 & sum(pmask=="sigma2")==0 ){ likelihoodOK = 0 } # theta_star[["error"]]<0 | theta_star[["tau1"]]<0 | theta_star[["tau2"]]<0 theta_star[["mu"]]<0 | theta_star[["mu"]]>10 |theta_star[["muShort"]]<0 | theta_star[["sigma"]]<0 | theta_star[["sigma2"]]<0  | theta_star[["wane"]]<0 | 
#   if(theta_star[["tau2"]]>10 & sum(pmask=="tau2")==0 ){ likelihoodOK = 0 } # theta_star[["error"]]<0 | theta_star[["tau1"]]<0 | theta_star[["tau2"]]<0 theta_star[["mu"]]<0 | theta_star[["mu"]]>10 |theta_star[["muShort"]]<0 | theta_star[["sigma"]]<0 | theta_star[["sigma2"]]<0  | theta_star[["wane"]]<0 | 


function AbstractMCMC.step(
    rng::Random.AbstractRNG,
    model::AbstractMCMC.LogDensityModel,
    sampler::MHParameterSampler,
    state::ParameterSamplerState;
    kwargs...
)
    prev_transition = state.transition
    theta = prev_transition.θ

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
    f = model.logdensity.ℓ

    context = DynamicPPL.DefaultContext()

    varinfo_prev = DynamicPPL.unflatten(f.varinfo, theta)
    logprob_previous = DynamicPPL.getlogp(last(DynamicPPL.evaluate!!(f.model, varinfo_prev, context)))

    # Theta proposal
    log_theta = log.(theta)

    theta_new_dist = MvNormal(log_theta, sigma_covar_adapted) # TODO add covar
    theta_new = exp.(rand(rng, theta_new_dist))

    # TODO add check if theta_new is out of bounds (per Kucharski model)

    # TODO replace with direct function call?
    varinfo_proposal = DynamicPPL.unflatten(f.varinfo, theta_new)
    logprob_proposal = DynamicPPL.getlogp(last(DynamicPPL.evaluate!!(f.model, varinfo_proposal, context)))

    if -Random.randexp(rng) <= logprob_proposal - logprob_previous
        # Accept theta_new
        transition = ParameterSamplerTransition(theta_new)

        n_accepted = n_accepted + 1
        return transition, ParameterSamplerState(transition, theta_new, sigma_covar_adapted, n_accepted, n_rejected)
    else
        # Reject theta_new, return theta
        transition = ParameterSamplerTransition(theta)

        n_rejected = n_rejected + 1
        return transition, ParameterSamplerState(transition, theta, sigma_covar_adapted, n_accepted, n_rejected)
    end
end

function make_mh_parameter_sampler()
    return externalsampler(MHParameterSampler(), adtype=Turing.DEFAULT_ADTYPE, unconstrained=false)
end