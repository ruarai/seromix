

function sample_chain(
    model,
    initial_params,
    gibbs_sampler,
    p,
    rng;
    n_sample::Int,
    n_thinning::Int,
    n_chain::Int,
    progress=true
)
    chain = sample(
        rng,
        model, gibbs_sampler, 
        MCMCThreads(), n_sample ÷ n_thinning, n_chain,

        thinning = n_thinning,
        callback = log_callback,
        initial_params = initial_params,
        progress = progress
    )

    check_inf_prob_birth_year(chain_infections_prob(chain, p), p)

    return chain
end

function model_symbols_apart_from(model, syms)
    symbols = DynamicPPL.syms(DynamicPPL.VarInfo(model))
    for s in syms
        symbols = symbols[findall(symbols .!= s)]
    end
    
    return symbols
end


function make_gibbs_sampler(model, p, proposal_function)
    symbols_not_inf = model_symbols_apart_from(model, [:infections])
    
    gibbs_sampler = Gibbs(
        :infections => make_mh_infection_sampler(p, proposal_function),
        symbols_not_inf => make_mh_parameter_sampler()
    )
    
    return gibbs_sampler
end


function log_callback(rng, model, sampler, sample, state, iteration; kwargs...)
    if iteration % 50 == 0
        if length(sampler.alg.samplers) > 1

            sampler_state = state.states[1].state

            rate = sampler_state.n_accepted / (sampler_state.n_accepted + sampler_state.n_rejected)

            time_elapsed = (sampler_state.time_B - sampler_state.time_A) * 1000

            # Time per 10,000 samples
            sampling_rate = time_elapsed / 6

            @printf("%6d; %6.2f; %7.2fms/sample;%7.2fmin/10,000 samples\n", iteration, rate, time_elapsed, sampling_rate)
        else
            @printf("%5d;\n", iteration)
        end
    end
end

function get_logp(theta, model; context = DynamicPPL.DefaultContext())
    f = model.logdensity
    varinfo_prev = DynamicPPL.unflatten(f.varinfo, theta)
    return DynamicPPL.getlogp(last(DynamicPPL.evaluate!!(f.model, varinfo_prev, context)))
end


# Hack to make sure logprob is carried through the Gibbs sampler
# function Turing.Inference.varinfo(state::Turing.Inference.TuringState)
#     if state.state.transition isa InfectionSamplerState
#         θ = state.state.transition.θ
#         vi = DynamicPPL.unflatten(state.ldf.varinfo, θ)
#         vi = setlogp!!(vi, state.state.transition.lp)
#     else
#         θ = Turing.Inference.getparams(state.ldf.model, state.state)
#         vi = DynamicPPL.unflatten(state.ldf.varinfo, θ)
#         return vi
#     end
# end

# Hack to reduce model evaluations
function Turing.Inference.transition_to_turing(f::DynamicPPL.LogDensityFunction, transition)
    if transition isa AdvancedMH.Transition
        θ = Turing.Inference.getparams(f.model, transition)
        varinfo = DynamicPPL.unflatten(f.varinfo, θ)
        return Turing.Inference.Transition(f.model, varinfo, transition)
    end



    θ = if transition isa InfectionSamplerTransition
        transition.θ
    elseif transition isa ParameterSamplerTransition
        transition.θ
    elseif transition isa SliceSampling.Transition
        Turing.Inference.getparams(f.model, transition) # TODO this may not work
    end

    varinfo = DynamicPPL.unflatten(f.varinfo, θ)
    return Turing.Inference.Transition(θ, transition.lp, Turing.Inference.getstats(transition))
end