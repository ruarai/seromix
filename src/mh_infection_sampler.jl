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

    # if haskey(kwargs, :initial_params)
    #     println("Setting initial infections matrix to initial_params")
    #     theta_init = kwargs[:initial_params]
    # else
    #     println("Setting initial infections matrix to false")
    #     theta_init = fill(false, d)
    # end

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
    prev_transition = state.transition
    theta = prev_transition.θ

    theta_new = copy(theta)

    # Get the logdensity function
    # println(fieldnames(typeof(model.logdensity)))
    f = model.logdensity

    varinfo_prev = DynamicPPL.unflatten(f.varinfo, theta)
    
    n_t_steps = sampler.n_t_steps
    n_subjects = sampler.n_subjects

    p_swap = 1.5 / n_t_steps

    mask = zeros(Bool, n_t_steps)
    # context = DynamicPPL.DefaultContext() ## TODO think about?

    for ix_subject in 1:n_subjects

        if rand(rng) < 0.8
            continue # Per Kucharski model, some chance of doing nothing for each individual
        end

        context = IndividualSubsetContext(ix_subject)
        logprob_previous = DynamicPPL.getlogp(last(DynamicPPL.evaluate!!(f.model, varinfo_prev, context)))

        # log_hastings_ratio = propose_mask_random!(rng, theta_new, mask, ix_subject, n_t_steps, p_swap)
        log_hastings_ratio = propose_mask_kucharski_literal!(rng, theta_new, mask, ix_subject, n_t_steps, p_swap)

        apply_mask!(theta_new, mask, ix_subject, n_t_steps)

        varinfo_proposal = DynamicPPL.unflatten(f.varinfo, theta_new)
        logprob_proposal = DynamicPPL.getlogp(last(DynamicPPL.evaluate!!(f.model, varinfo_proposal, context)))

        if -Random.randexp(rng) <= logprob_proposal - logprob_previous + log_hastings_ratio
            # Accept, do nothing

            sampler.acceptions += 1
        else
            # Reject, re-apply mask to undo
            apply_mask!(theta_new, mask, ix_subject, n_t_steps)
            sampler.rejections += 1
        end
    end


    transition = InfectionSamplerTransition(theta_new)
    return transition, InfectionSamplerState(transition, theta_new)
end

function propose_mask_random!(
    rng, theta::AbstractVector{Bool}, mask::Vector{Bool}, ix_subject::Int, n_t_steps::Int,
    p_swap::Real
)
    mask .= rand(rng, Bernoulli(p_swap), n_t_steps)

    return 0
end

function propose_mask_kucharski_literal!(
    rng, theta::AbstractVector{Bool}, mask::Vector{Bool}, ix_subject::Int, n_t_steps::Int,
    p_swap::Real
)
    mask .= false

    ix_start = (ix_subject - 1) * n_t_steps + 1
    ix_end = ix_start + n_t_steps - 1

    r_sample = rand(rng)

    if r_sample < 1/3
        # Case 1 - remove an infection

        inf_indices = findall(@view theta[ix_start:ix_end])

        if length(inf_indices) > 0
            ix_remove = sample(rng, inf_indices)
            mask[ix_remove] = true

            # Forward transition probability
            # 1/3 * 1 / n_I(s)
            p_forward = 1/3 * 1 / length(inf_indices)

            # Reverse transition probability, from case 2
            # 1/3 * 1 / n_U(s+1) = 1/3 * 1 / (n_t_steps - n_I(s) + 1)
            p_reverse = 1/3 * 1 / (n_t_steps - length(inf_indices) + 1)

            return log(p_reverse) - log(p_forward) # TODO simplify out the log terms?
        end
    elseif r_sample < 2/3
        # Case 2 - add an infection

        not_inf_indices = findall(.!(@view theta[ix_start:ix_end]))

        if length(not_inf_indices) > 0
            ix_remove = sample(rng, not_inf_indices)
            mask[ix_remove] = true

            # Forward transition probability
            # 1/3 * 1 / n_U(s)
            p_forward = 1/3 * 1 / length(not_inf_indices)

            # Reverse transition probability, from case 1
            # 1/3 * 1 / n_I(s+1) = 1/3 * 1 / (n_t_steps - n_U(s) + 1)
            p_reverse = 1/3 * 1 / (n_t_steps - length(not_inf_indices) + 1)

            return log(p_reverse) - log(p_forward)
        end
    else
        # Case 3 - move an infection

        inf_indices = findall(@view theta[ix_start:ix_end])
        not_inf_indices = findall(.!(@view theta[ix_start:ix_end]))

        if length(inf_indices) > 0 && length(not_inf_indices) > 0
            ix_t_from = sample(rng, inf_indices)
            ix_t_to = sample(rng, not_inf_indices)

            mask[ix_t_from] = true
            mask[ix_t_to] = true

            # Forward and reverse transition probabilities are equal
            return 0
        end
    end

    return 0 # Nothing has occurred, no hastings ratio
end

function apply_mask!(theta::AbstractVector{Bool}, mask::Vector{Bool}, ix_subject::Int, n_t_steps::Int)
    ix_start = (ix_subject - 1) * n_t_steps + 1
    ix_end = ix_start + n_t_steps - 1
    
    @inbounds for (i, i_theta) in enumerate(ix_start:ix_end)
        theta[i_theta] = xor(theta[i_theta], mask[i])
    end
end

struct IndividualSubsetContext <: DynamicPPL.AbstractContext
    ix::Int
end

DynamicPPL.NodeTrait(context::IndividualSubsetContext) = DynamicPPL.IsLeaf()

function make_mh_infection_sampler(n_t_steps, n_subjects)
    return externalsampler(MHInfectionSampler(n_t_steps, n_subjects), adtype=Turing.DEFAULT_ADTYPE, unconstrained=false)
end