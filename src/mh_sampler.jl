mutable struct MHInfectionSampler <: AbstractMCMC.AbstractSampler
end

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
    spl::MHInfectionSampler;
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
    
    # TODO --- how to get n_subjects here?
    # and n_t_steps?
    n_t_steps = 20
    n_subjects = 20

    mean_p_swap = 1.0 / n_t_steps
    alpha = 3.0

    p_swap = rand(rng, Beta(alpha, alpha / mean_p_swap - alpha))

    mask = zeros(Bool, n_t_steps)

    ## TODO --- verify this individual-by-individual sampling
    # is mathematically correct
    for ix_subject in 1:n_subjects
        context = IndividualSubsetContext(ix_subject)
        
        logprob_previous = DynamicPPL.getlogp(last(DynamicPPL.evaluate!!(f.model, varinfo_prev, context)))

        mask .= rand(rng, Bernoulli(p_swap), n_t_steps)

        apply_mask!(theta_new, mask, ix_subject, n_t_steps)

        varinfo_proposal = DynamicPPL.unflatten(f.varinfo, theta_new)

        logprob_proposal = DynamicPPL.getlogp(last(DynamicPPL.evaluate!!(f.model, varinfo_proposal, context)))

        if -Random.randexp(rng) <= logprob_proposal - logprob_previous
            # nothing
        else
            # Re-apply mask to undo
            apply_mask!(theta_new, mask, ix_subject, n_t_steps)
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

function make_mh_infection_sampler()
    return externalsampler(MHInfectionSampler(), adtype=Turing.DEFAULT_ADTYPE, unconstrained=false)
end


function DynamicPPL.unflatten(vi::DynamicPPL.VarInfo, spl::AbstractMCMC.AbstractSampler, x::AbstractVector)
    md = DynamicPPL.unflatten(vi.metadata, spl, x)

    return DynamicPPL.VarInfo(
        md,
        Base.RefValue{DynamicPPL.float_type_with_fallback(eltype(x))}(DynamicPPL.getlogp(vi)),
        Ref(DynamicPPL.get_num_produce(vi)),
    )
end

# function DynamicPPL.unflatten(vi::DynamicPPL.VarInfo, x::AbstractVector)
#     md = DynamicPPL.unflatten(vi.metadata, spl, x)
#     # Note that use of RefValue{eltype(x)} rather than Ref is necessary to deal with cases
#     # where e.g. x is a type gradient of some AD backend.
#     return DynamicPPL.VarInfo(
#         md,
#         Base.RefValue{float_type_with_fallback(eltype(x))}(DynamicPPL.getlogp(vi)),
#         Ref(DynamicPPL.get_num_produce(vi)),
#     )
# end
