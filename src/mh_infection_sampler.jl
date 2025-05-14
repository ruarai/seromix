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

    theta_init = fill(false, d)

    # theta_init = rand(rng, Bernoulli(0.1), d)

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
    f = model.logdensity.ℓ

    varinfo_prev = DynamicPPL.unflatten(f.varinfo, theta)
    
    n_t_steps = sampler.n_t_steps
    n_subjects = sampler.n_subjects

    p_swap = 1.5 / n_t_steps

    mask = zeros(Bool, n_t_steps)
    # context = DynamicPPL.DefaultContext() ## TODO think about?

    for ix_subject in 1:n_subjects
        context = IndividualSubsetContext(ix_subject)
        logprob_previous = DynamicPPL.getlogp(last(DynamicPPL.evaluate!!(f.model, varinfo_prev, context)))

        propose_mask_random!(rng, theta_new, mask, ix_subject, n_t_steps, p_swap)

        apply_mask!(theta_new, mask, ix_subject, n_t_steps)

        varinfo_proposal = DynamicPPL.unflatten(f.varinfo, theta_new)
        logprob_proposal = DynamicPPL.getlogp(last(DynamicPPL.evaluate!!(f.model, varinfo_proposal, context)))

        if -Random.randexp(rng) <= logprob_proposal - logprob_previous
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

function propose_mask_random!(rng, theta::AbstractVector{Bool}, mask::Vector{Bool}, ix_subject::Int, n_t_steps::Int, p_swap::Real)

    mask .= rand(rng, Bernoulli(p_swap), n_t_steps)

end

function propose_mask_kucharski!(rng, theta::AbstractVector{Bool}, mask::Vector{Bool}, ix_subject::Int, n_t_steps::Int, p_swap::Real)

    mask .= false

    if rand(rng) < 0.5
        return # 50% chance of doing nothing
    end

    ix_start = (ix_subject - 1) * n_t_steps + 1
    ix_end = ix_start + n_t_steps - 1

    if rand(rng) < 0.66
        # Case 1 - remove/add an infection

        ix_swap = sample(rng, 1:n_t_steps)
        mask[ix_swap] = true
    else
        # Case 3 - move an infection

        n_inf = sum(@view theta[ix_start:ix_end])

        if n_inf > 0
            # Maybe just sample randomly and check?
            inf_indices = findall(@view theta[ix_start:ix_end])

            ix_t_from = sample(rng, inf_indices)
            ix_t_to = sample(rng, 1:n_t_steps)

            # Ensure ix_t_to is not currently an infection event
            if !view(theta, ix_start:ix_end)[ix_t_to]
                # Set mask such that ix_t_from will become not an infection event
                # and ix_t_to will be an infection event
                # i.e. move infection
                mask[ix_t_from] = true
                mask[ix_t_to] = true
            end
        end
    end
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