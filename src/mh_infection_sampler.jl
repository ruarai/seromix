mutable struct MHInfectionSampler <: AbstractMCMC.AbstractSampler
    # Ideally these would just be in the model but not sure how to do that.
    n_t_steps::Int
    n_subjects::Int
    
    # TODO move these out of here
    rejections::Int
    acceptions::Int
end

MHInfectionSampler(n_t_step, n_subjects) = MHInfectionSampler(n_t_step, n_subjects, 0, 0)

mh_sampler_acceptance_rate(s::MHInfectionSampler) = s.acceptions / (s.acceptions + s.rejections)

isgibbscomponent(::MHInfectionSampler) = true

struct InfectionSamplerTransition{T}
    θ::AbstractVector{T}
end

struct InfectionSamplerState{P<:InfectionSamplerTransition, T}
    transition::P
    θ::AbstractVector{T}
end

AbstractMCMC.getparams(state::InfectionSamplerState) = state.θ
function AbstractMCMC.setparams!!(state::InfectionSamplerState, θ)
    return InfectionSamplerState(InfectionSamplerTransition(θ), θ)
end


function AbstractMCMC.step(
    rng::Random.AbstractRNG,
    model::AbstractMCMC.LogDensityModel,
    sampler::MHInfectionSampler;
    kwargs...
)
    d = LogDensityProblems.dimension(model.logdensity)

    theta_init = rand(rng, Bernoulli(0.1), d)

    transition = InfectionSamplerTransition(theta_init)
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

    for ix_subject in 1:n_subjects

        if rand(rng) < 0.8
            continue # Per Kucharski model, some chance of doing nothing for each individual
        end

        # Set the context of DynamicPPL to indicate that we are only
        # calculating likelihood over one individual. Note will
        # still calculate prior across the infection history matrix
        context = IndividualSubsetContext(ix_subject)
        logprob_previous = DynamicPPL.getlogp(last(DynamicPPL.evaluate!!(f.model, varinfo_prev, context)))

        swap_indices, log_hastings_ratio = propose_swaps_v2!(rng, theta_new, ix_subject, n_t_steps)

        apply_swaps!(theta_new, swap_indices, ix_subject, n_t_steps)

        varinfo_proposal = DynamicPPL.unflatten(f.varinfo, theta_new)
        logprob_proposal = DynamicPPL.getlogp(last(DynamicPPL.evaluate!!(f.model, varinfo_proposal, context)))

        if -Random.randexp(rng) <= logprob_proposal - logprob_previous + log_hastings_ratio
            # Accept, do nothing

            sampler.acceptions += 1
        else
            # Reject, re-apply swaps to undo
            apply_swaps!(theta_new, swap_indices, ix_subject, n_t_steps)
            sampler.rejections += 1
        end
    end


    transition = InfectionSamplerTransition(theta_new)
    return transition, InfectionSamplerState(transition, theta_new)
end

function apply_swaps!(theta, swap_indices, ix_subject, n_t_steps)
    ix_start = (ix_subject - 1) * n_t_steps

    for ix_swap in swap_indices
        theta[ix_start + ix_swap] = !theta[ix_start + ix_swap]
    end
end


struct IndividualSubsetContext <: DynamicPPL.AbstractContext
    ix::Int
end

DynamicPPL.NodeTrait(context::IndividualSubsetContext) = DynamicPPL.IsLeaf()

function make_mh_infection_sampler(n_t_steps, n_subjects)
    return externalsampler(MHInfectionSampler(n_t_steps, n_subjects), adtype=Turing.DEFAULT_ADTYPE, unconstrained=false)
end