

# A re-implementation of the infection history sampler as
# used in Kucharski (2018)

# For each individual, samples a new history using step_function

mutable struct MHInfectionSampler <: AbstractMCMC.AbstractSampler
    # Ideally these would just be in the model but not sure how to do that.
    n_t_steps::Int
    n_subjects::Int

    prop_sample::Float64

    step_function
end

struct InfectionSamplerTransition{T, L <: Real}
    θ::AbstractVector{T}
    lp::L
end

struct InfectionSamplerState{P<:InfectionSamplerTransition, T}
    transition::P
    θ::AbstractVector{T}

    n_accepted::Int
    n_rejected::Int

    time_A::Float64
    time_B::Float64
end

MHInfectionSampler(n_t_step, n_subjects, step_fn) = MHInfectionSampler(n_t_step, n_subjects, 1.0, step_fn)

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
    return transition, InfectionSamplerState(transition, theta_init, 0, 0, time(), time())
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
    # subject_indices = sample(rng, 1:n_subjects, ceil(Int, 0.4 * n_subjects))
    # TODO add as option.

    subject_indices = 1:n_subjects

    accepted = 0
    rejected = 0

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

            accepted += 1
        else
            # Reject, re-apply swaps to undo

            rejected += 1

            apply_swaps!(theta_new, swap_indices, ix_subject, n_t_steps)
        end
    end


    transition = InfectionSamplerTransition(theta_new, -Inf)

    return transition, InfectionSamplerState(transition, theta_new, accepted, rejected, state.time_B, time())
end

function apply_swaps!(theta, swap_indices, ix_subject, n_t_steps)
    ix_start = (ix_subject - 1) * n_t_steps

    @inbounds for ix_swap in swap_indices
        theta[ix_start + ix_swap] = !theta[ix_start + ix_swap]
    end
end


# Below is necessary framework for use in Turing

# Puts the sampler into an external sampler for use in Turing
function make_mh_infection_sampler(n_t_steps, n_subjects, step_fn)
    return externalsampler(MHInfectionSampler(n_t_steps, n_subjects, step_fn), adtype=Turing.DEFAULT_ADTYPE, unconstrained=false)
end

isgibbscomponent(::MHInfectionSampler) = true

AbstractMCMC.getparams(state::InfectionSamplerState) = state.θ
function AbstractMCMC.setparams!!(state::InfectionSamplerState, θ)
    return InfectionSamplerState(InfectionSamplerTransition(θ, state.transition.lp), θ, state.n_accepted, state.n_rejected, state.time_A, state.time_B)
end


# A DynamicPPL context that indicates to the model that we are only
# calculating likelihood over one individual (subject)
struct IndividualSubsetContext <: DynamicPPL.AbstractContext
    ix::Int
end

DynamicPPL.NodeTrait(context::IndividualSubsetContext) = DynamicPPL.IsLeaf()
