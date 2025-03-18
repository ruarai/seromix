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

    theta_init = convert(Vector{Real}, fill(false, d))

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

    # TOOD -- optimise the code here. lots of copying at the moment

    # TODO --- how to get n_ind here?
    # and n_t_steps?

    ## TODO --- verify this individual-by-individual sampling
    # is mathematically correct
    for ix_ind in 1:69
        context = IndividualSubsetContext(ix_ind)
        
        logprob_previous = DynamicPPL.getlogp(last(DynamicPPL.evaluate!!(f.model, varinfo_prev, context)))

        theta_proposal = propose_theta(rng, theta_new, ix_ind)

        varinfo_proposal = DynamicPPL.unflatten(f.varinfo, theta_proposal)

        logprob_proposal = DynamicPPL.getlogp(last(DynamicPPL.evaluate!!(f.model, varinfo_proposal, context)))

        if -Random.randexp(rng) <= logprob_proposal - logprob_previous
            theta_new = theta_proposal
        else
            # nothing?
        end
    end


    transition = Transition(theta_new)
    return transition, SamplerState(transition, theta_new)
end

function propose_theta(rng, theta::AbstractVector{T}, ix_ind::Int) where T <: Real
    n_t_steps = 45
    p_swap = rand(rng, Beta(2, n_t_steps - 2))
    
    # Generate mask directly (single allocation)
    theta_mask = rand(rng, Bernoulli(p_swap), n_t_steps)
    theta_prop = copy(theta)

    ix_start = (ix_ind - 1) * n_t_steps + 1
    ix_end = ix_start + n_t_steps - 1
    
    @inbounds for (i, i_theta) in enumerate(ix_start:ix_end)
        theta_prop[i_theta] = xor(theta[i_theta], theta_mask[i])
    end
    
    return Vector{Real}(theta_prop)
end


struct IndividualSubsetContext <: DynamicPPL.AbstractContext
    ix::Int
end

DynamicPPL.NodeTrait(context::IndividualSubsetContext) = DynamicPPL.IsLeaf()

function make_mh_infection_sampler()
    return externalsampler(MHInfectionSampler(), adtype=Turing.DEFAULT_ADTYPE, unconstrained=false)
end