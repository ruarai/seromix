

function sample_chain(
    model, 
    gibbs_sampler;
    n_sample::Int,
    n_thinning::Int,
    n_chain::Int
)
    return sample(
        model, gibbs_sampler, 
        MCMCThreads(), n_sample ÷ n_thinning, n_chain,

        thinning = n_thinning,
        callback = log_callback
    )
end

function model_symbols_apart_from(model, sym)
    symbols = DynamicPPL.syms(DynamicPPL.VarInfo(model))
    symbols = symbols[findall(symbols .!= sym)]
    
    return symbols
end

function make_gibbs_sampler(model, inf_sym, hmc_step_size, n_leapfrog, params)
    symbols_not_inf = model_symbols_apart_from(model, inf_sym)
    
    # Must somehow balance the level of exploration of the MH sampler
    # with that of the HMC sampler -- so repeating MH or changing HMC step size
    gibbs_sampler = Gibbs(
        :infections => make_mh_infection_sampler(params.n_t_steps, params.n_subjects),
        symbols_not_inf => HMC(hmc_step_size, n_leapfrog) # Must be reduced with number of individuals?
    )
    
    return gibbs_sampler
end


function make_gibbs_sampler(model, inf_sym)
    symbols_not_inf = model_symbols_apart_from(model, inf_sym)
    
    gibbs_sampler = Gibbs(
        :infections => make_mh_infection_sampler(p.n_t_steps, p.n_subjects),
        symbols_not_inf => make_mh_parameter_sampler()
    )
    
    return gibbs_sampler
end

function log_callback(rng, model, sampler, sample, state, iteration; kwargs...)
    
    if iteration % 50 == 0
        # print("$iteration,")
        if length(sampler.alg.samplers) > 1
            mh_sampler = sampler.alg.samplers[1].alg.sampler

            # TODO - get HMC acceptance rate here?
            # Is the current method weird, also?

            println("$iteration, $(mh_sampler_acceptance_rate(mh_sampler))")

            mh_sampler.acceptions = 0
            mh_sampler.rejections = 0
        else
            println("$iteration")
        end
    end
end