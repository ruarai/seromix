

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
            mh_sampler = sampler.alg.samplers[1].alg.sampler

            rate = mh_sampler_acceptance_rate(mh_sampler)
            println("$iteration, $(round(rate, digits=2))")

            mh_sampler.acceptions = 0
            mh_sampler.rejections = 0
        else
            println("$iteration")
        end
    end
end


# Hack to make sure logprob is carried through the Gibbs sampler
function Turing.Inference.varinfo(state::Turing.Inference.TuringState)
    θ = Turing.Inference.getparams(state.ldf.model, state.state)
    vi = DynamicPPL.unflatten(state.ldf.varinfo, θ)

    vi = setlogp!!(vi, state.state.transition.lp)

    return vi
end
