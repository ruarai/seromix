

function sample_chain(
    model,
    initial_params,
    gibbs_sampler,
    rng;
    n_sample::Int,
    n_thinning::Int,
    n_chain::Int,
    progress=true
)
    return sample(
        rng,
        model, gibbs_sampler, 
        MCMCThreads(), n_sample ÷ n_thinning, n_chain,

        thinning = n_thinning,
        callback = log_callback,
        initial_params = initial_params,
        progress = progress
    )
end

function model_symbols_apart_from(model, syms)
    symbols = DynamicPPL.syms(DynamicPPL.VarInfo(model))
    for s in syms
        symbols = symbols[findall(symbols .!= s)]
    end
    
    return symbols
end


function make_gibbs_sampler(model, p, step_fn)
    symbols_not_inf = model_symbols_apart_from(model, [:infections])
    
    gibbs_sampler = Gibbs(
        :infections => make_mh_infection_sampler(p.n_t_steps, p.n_subjects, step_fn),
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

            time_for_thousand =  time_elapsed / 60

            # println("$iteration; $(round(rate, digits=2)); $(round(time_elapsed, digits=2))ms.")

            @printf("%6d; %6.2f; %7.2fms;%7.2fmin/1,000\n", iteration, rate, time_elapsed, time_for_thousand)
        else
            @printf("%5d;\n", iteration)
        end
    end
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
# will likely break anything with other samplers, bijection
# function Turing.Inference.transition_to_turing(f::DynamicPPL.LogDensityFunction, transition)
#     θ = transition.θ
#     varinfo = DynamicPPL.unflatten(f.varinfo, θ)
#     return Turing.Inference.Transition(θ, transition.lp, Turing.Inference.getstats(transition))
# end