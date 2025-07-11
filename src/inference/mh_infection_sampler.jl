

# A re-implementation of the infection history sampler as
# used in Kucharski (2018)

# For each individual, samples a new history using proposal_function

mutable struct MHInfectionSampler2 <: AbstractMCMC.AbstractSampler
    n_t_steps::Int
    n_subjects::Int

    subject_birth_ix::Vector{Int}

    proposal_proportion::Float64
    n_repeats::Int

    proposal_function::Function
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

MHInfectionSampler2(n_t_step, n_subjects, subject_birth_ix, proposal_function) = MHInfectionSampler2(n_t_step, n_subjects, subject_birth_ix, 1.0, proposal_function)

function AbstractMCMC.step(
    rng::Random.AbstractRNG,
    model::AbstractMCMC.LogDensityModel,
    sampler::MHInfectionSampler2;
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

    logprob_init = get_logp(theta_init, model)

    transition = InfectionSamplerTransition(theta_init, logprob_init)
    return transition, InfectionSamplerState(transition, theta_init, 0, 0, time(), time())
end

function AbstractMCMC.step(
    rng::Random.AbstractRNG,
    model::AbstractMCMC.LogDensityModel,
    sampler::MHInfectionSampler2,
    state::InfectionSamplerState;
    kwargs...
)
    theta = state.transition.θ
    theta_new = copy(theta) # TODO check if this can be removed

    n_t_steps = sampler.n_t_steps
    n_subjects = sampler.n_subjects

    accepted = 0
    rejected = 0

    # Using ceil rather than round here to match original model precisely
    n_sample = ceil(Int, sampler.proposal_proportion * n_subjects)
    # Buffer to hold subject indices
    subject_indices = Vector{Int}(undef, n_sample)
    
    # If n_sample is equal to n_subjects, then we work across 1:n_subjects
    if n_sample == n_subjects
        subject_indices .= 1:n_subjects
    end

    for _ in 1:sampler.n_repeats
        # Fill subject_indices in-place with random sample (without replacement)
        if n_sample != n_subjects
            sample!(rng, 1:n_subjects, subject_indices, replace = false)
        end

        for ix_subject in subject_indices
            # Set the context of DynamicPPL to indicate that we are only
            # calculating likelihood over one individual. Note will
            # still calculate prior across the infection history matrix
            context = IndividualSubsetContext(ix_subject)
            logprob_previous = get_logp(theta_new, model; context = context)

            subject_birth_ix = sampler.subject_birth_ix[ix_subject]

            swap_indices, log_hastings_ratio = sampler.proposal_function(rng, theta_new, ix_subject, n_t_steps, subject_birth_ix)

            apply_swaps!(theta_new, swap_indices, ix_subject, n_t_steps, subject_birth_ix)

            logprob_proposal = get_logp(theta_new, model; context = context)

            if -Random.randexp(rng) <= logprob_proposal - logprob_previous + log_hastings_ratio
                # Accept, do nothing

                accepted += 1
            else
                # Reject, re-apply swaps to undo

                rejected += 1
                apply_swaps!(theta_new, swap_indices, ix_subject, n_t_steps, subject_birth_ix)
            end
        end
    end

    transition = InfectionSamplerTransition(theta_new, -Inf)

    return transition, InfectionSamplerState(transition, theta_new, accepted, rejected, state.time_B, time())
end

function apply_swaps!(theta, swap_indices, ix_subject, n_t_steps, subject_birth_ix)
    ix_start = (ix_subject - 1) * n_t_steps + max(1, subject_birth_ix) - 1

    @inbounds for ix_swap in swap_indices
        theta[ix_start + ix_swap] = !theta[ix_start + ix_swap]
    end
end


# Below is necessary framework for use in Turing

# Puts the sampler into an external sampler for use in Turing
function make_mh_infection_sampler(sp, proposal_function; prop_sample = 0.4, n_repeats = 1)
    return externalsampler(
        MHInfectionSampler2(
            sp.n_t_steps,sp.n_subjects, 
            sp.subject_birth_ix,
            prop_sample,
            n_repeats,
            proposal_function
        ), 
        adtype=Turing.DEFAULT_ADTYPE,
        unconstrained=false
    )
end


isgibbscomponent(::MHInfectionSampler2) = true

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
