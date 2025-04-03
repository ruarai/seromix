mutable struct MHInfectionSampler <: AbstractMCMC.AbstractSampler
    # Ideally these would just be in the model but not sure how to do that.
    n_t_steps::Int
    n_subjects::Int
    
    rejections::Int
    acceptions::Int
end

MHInfectionSampler(n_t_step, n_subjects) = MHInfectionSampler(n_t_step, n_subjects, 0, 0)

mh_sampler_acceptance_rate(s::MHInfectionSampler) = s.acceptions / (s.acceptions + s.rejections)

import Turing.Inference: isgibbscomponent

isgibbscomponent(::MHInfectionSampler) = true

struct Transition{T}
    θ::AbstractVector{T}
end

struct SamplerState{P<:Transition, T}
    transition::P
    θ::AbstractVector{T}
end

AbstractMCMC.getparams(state::SamplerState) = state.θ
function AbstractMCMC.setparams!!(state::SamplerState, θ)
    return SamplerState(
        Transition(θ), θ
    )
end


function AbstractMCMC.step(
    rng::Random.AbstractRNG,
    model::AbstractMCMC.LogDensityModel,
    sampler::MHInfectionSampler;
    kwargs...
)
    d = LogDensityProblems.dimension(model.logdensity)

    # theta_init = convert(Vector{Real}, fill(false, d))
    theta_init = fill(false, d)

    transition = Transition(theta_init)
    return transition, SamplerState(transition, theta_init)
end

function AbstractMCMC.step(
    rng::Random.AbstractRNG,
    model::AbstractMCMC.LogDensityModel,
    sampler::MHInfectionSampler,
    state::SamplerState;
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

    p_swap = 0.1 / n_t_steps

    mask = zeros(Bool, n_t_steps)

    for ix_subject in 1:n_subjects
        context = IndividualSubsetContext(ix_subject)
        
        logprob_previous = DynamicPPL.getlogp(last(DynamicPPL.evaluate!!(f.model, varinfo_prev, context)))

        mask .= rand(rng, Bernoulli(p_swap), n_t_steps)

        apply_mask!(theta_new, mask, ix_subject, n_t_steps)

        varinfo_proposal = DynamicPPL.unflatten(f.varinfo, theta_new)

        logprob_proposal = DynamicPPL.getlogp(last(DynamicPPL.evaluate!!(f.model, varinfo_proposal, context)))

        if -Random.randexp(rng) <= logprob_proposal - logprob_previous
            # nothing

            sampler.acceptions += 1
        else
            # Re-apply mask to undo
            apply_mask!(theta_new, mask, ix_subject, n_t_steps)
            sampler.rejections += 1
        end
    end


    transition = Transition(theta_new)
    return transition, SamplerState(transition, theta_new)
end

function apply_mask!(theta::AbstractVector{Bool}, mask::Vector{Bool}, ix_ind::Int, n_t_steps::Int)
    ix_start = (ix_ind - 1) * n_t_steps + 1
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


function DynamicPPL.unflatten(vi::DynamicPPL.VarInfo, spl::AbstractMCMC.AbstractSampler, x::AbstractVector)
    md = DynamicPPL.unflatten(vi.metadata, spl, x)

    return DynamicPPL.VarInfo(
        md,
        Base.RefValue{DynamicPPL.float_type_with_fallback(eltype(x))}(DynamicPPL.getlogp(vi)),
        Ref(DynamicPPL.get_num_produce(vi)),
    )
end