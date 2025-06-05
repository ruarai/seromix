mutable struct MHInfectionSampler <: AbstractMCMC.AbstractSampler
    # Ideally these would just be in the model but not sure how to do that.
    n_t_steps::Int
    n_subjects::Int
    
    # TODO move these out of here
    rejections::Int
    acceptions::Int

    step_function
end

MHInfectionSampler(n_t_step, n_subjects, step_fn) = MHInfectionSampler(n_t_step, n_subjects, 0, 0, step_fn)

mh_sampler_acceptance_rate(s::MHInfectionSampler) = s.acceptions / (s.acceptions + s.rejections)

isgibbscomponent(::MHInfectionSampler) = true

struct InfectionSamplerTransition{T, L <: Real}
    θ::AbstractVector{T}
    lp::L
end

struct InfectionSamplerState{P<:InfectionSamplerTransition, T}
    transition::P
    θ::AbstractVector{T}
end

AbstractMCMC.getparams(state::InfectionSamplerState) = state.θ
function AbstractMCMC.setparams!!(state::InfectionSamplerState, θ)
    return InfectionSamplerState(InfectionSamplerTransition(θ, state.transition.lp), θ)
end


function AbstractMCMC.step(
    rng::Random.AbstractRNG,
    model::AbstractMCMC.LogDensityModel,
    sampler::MHInfectionSampler;
    initial_params,
    kwargs...
)
    d = LogDensityProblems.dimension(model.logdensity)

    if !isnothing(initial_params)
        println("Setting initial infections matrix to initial_params")
        theta_init = initial_params
    else
        println("Setting initial infections matrix to false")
        theta_init = fill(false, d)
    end

    varinfo_init = DynamicPPL.unflatten(model.logdensity.varinfo, theta_init)
    logprob_init = DynamicPPL.getlogp(last(DynamicPPL.evaluate!!(model.logdensity.model, varinfo_init, DynamicPPL.DefaultContext())))

    transition = InfectionSamplerTransition(theta_init, logprob_init)
    return transition, InfectionSamplerState(transition, theta_init)
end

function AbstractMCMC.step(
    rng::Random.AbstractRNG,
    model::AbstractMCMC.LogDensityModel,
    sampler::MHInfectionSampler,
    state::InfectionSamplerState;
    kwargs...
)
    theta = state.transition.θ
    theta_new = copy(theta)

    # Get the logdensity function
    f = model.logdensity

    varinfo_prev = DynamicPPL.unflatten(f.varinfo, theta)
    
    n_t_steps = sampler.n_t_steps
    n_subjects = sampler.n_subjects

    # Per Kucharski model, only step for some % of individuals
    subject_indices = sample(rng, 1:n_subjects, floor(Int, 0.4 * n_subjects))

    logprob_final = 0.0

    for ix_subject in subject_indices
        # Set the context of DynamicPPL to indicate that we are only
        # calculating likelihood over one individual. Note will
        # still calculate prior across the infection history matrix
        context = IndividualSubsetContext(ix_subject)
        logprob_previous = DynamicPPL.getlogp(last(DynamicPPL.evaluate!!(f.model, varinfo_prev, context)))

        swap_indices, log_hastings_ratio = sampler.step_function(rng, theta_new, ix_subject, n_t_steps)

        apply_swaps!(theta_new, swap_indices, ix_subject, n_t_steps)

        varinfo_proposal = DynamicPPL.unflatten(f.varinfo, theta_new)
        logprob_proposal = DynamicPPL.getlogp(last(DynamicPPL.evaluate!!(f.model, varinfo_proposal, context)))

        if -Random.randexp(rng) <= logprob_proposal - logprob_previous + log_hastings_ratio
            # Accept, do nothing

            logprob_final = logprob_proposal
            sampler.acceptions += 1
        else
            # Reject, re-apply swaps to undo

            logprob_final = logprob_previous

            apply_swaps!(theta_new, swap_indices, ix_subject, n_t_steps)
            sampler.rejections += 1
        end
    end


    transition = InfectionSamplerTransition(theta_new, logprob_final)
    return transition, InfectionSamplerState(transition, theta_new)
end

function apply_swaps!(theta, swap_indices, ix_subject, n_t_steps)
    ix_start = (ix_subject - 1) * n_t_steps

    @inbounds for ix_swap in swap_indices
        theta[ix_start + ix_swap] = !theta[ix_start + ix_swap]
    end
end


struct IndividualSubsetContext <: DynamicPPL.AbstractContext
    ix::Int
end

DynamicPPL.NodeTrait(context::IndividualSubsetContext) = DynamicPPL.IsLeaf()

function make_mh_infection_sampler(n_t_steps, n_subjects, step_fn)
    return externalsampler(MHInfectionSampler(n_t_steps, n_subjects, step_fn), adtype=Turing.DEFAULT_ADTYPE, unconstrained=false)
end