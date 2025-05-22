

function sample_chain(
    model,
    initial_params,
    gibbs_sampler,
    rng;
    n_sample::Int,
    n_thinning::Int,
    n_chain::Int
)
    return sample(
        rng,
        model, gibbs_sampler, 
        MCMCThreads(), n_sample ÷ n_thinning, n_chain,

        thinning = n_thinning,
        callback = log_callback,
        initial_params = initial_params
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
    symbols_not_inf = model_symbols_apart_from(model, :infections)
    
    gibbs_sampler = Gibbs(
        :infections => make_mh_infection_sampler(p.n_t_steps, p.n_subjects, step_fn),
        symbols_not_inf => make_mh_parameter_sampler()
    )
    
    return gibbs_sampler
end

function make_initial_params_data_study(n_chain, init_matrix, rng)
    return [(
        mu_long = 2.0 + rand(rng, Uniform(-0.2, 0.2)),
        mu_short = 2.5 + rand(rng, Uniform(-0.2, 0.2)), 
        omega = 0.8 + rand(rng, Uniform(-0.05, 0.05)), 
        sigma_long = 0.15 + rand(rng, Uniform(-0.02, 0.02)),
        sigma_short = 0.05 + rand(rng, Uniform(-0.005, 0.005)), 
        tau = 0.05 + rand(rng, Uniform(-0.01, 0.01)), 
        obs_sd = 1.5 + rand(rng, Uniform(-0.1, 0.1)), 

        infections = init_matrix
    ) for i in 1:n_chain]
end

function make_initial_params_sim_study(p, obs_df, n_chain, rng)

    param_link_dists = [
        Uniform(0, 10), Uniform(0, 10),
        Uniform(0, 1), Uniform(0, 10),
        Uniform(0, 10), Uniform(0, 10),
        Uniform(0, 10)
    ]
    param_means = [2.0, 2.5, 0.8, 0.15, 0.05, 0.05, 1.0]


    return [(
        params = [link(param_link_dists[i], param_means[i]) for i in 1:7], 
        infections = initial_infections_matrix(p, obs_df, rng)
    ) for i in 1:n_chain]
end

# Translated from Kucharski model, but not sure that it's necessary.
function initial_infections_matrix(p, obs_df, rng)
    infections_0 = zeros(Bool, p.n_t_steps, p.n_subjects)

    obs_df_grouped = groupby(obs_df, :ix_subject)

    for ix_subject in 1:p.n_subjects
        
        obs_df_subject = obs_df_grouped[ix_subject]

        # Any observations where titre > 4 start as an infection
        for obs_row in eachrow(obs_df_subject)
            if obs_row.observed_titre > 4.0 && rand(rng) > 0.1
                infections_0[obs_row.ix_strain, ix_subject] = true
            end
        end

        # Make sure there is at least one infection in the 5 years preceding a test
        for obs_row in eachrow(obs_df_subject)
            t_range = max(1, obs_row.ix_t_obs - 5):obs_row.ix_t_obs
            if sum(infections_0[t_range, ix_subject]) == 0
                ix_t_sample = sample(rng, t_range)

                infections_0[ix_t_sample, ix_subject] = true
            end
        end

        # Make sure an individual only has a maximum of 10 infections (idk why)
        while sum(infections_0[:, ix_subject]) > 10
            ix_inf_drop = sample(rng, findall(infections_0[:, ix_subject]))
            infections_0[ix_inf_drop, ix_subject] = false
        end
    end

    mask_infections_birth_year!(infections_0, p.subject_birth_ix)

    return infections_0
end

function log_callback(rng, model, sampler, sample, state, iteration; kwargs...)
    
    if iteration % 50 == 0
        # print("$iteration,")
        if length(sampler.alg.samplers) > 1
            mh_sampler = sampler.alg.samplers[1].alg.sampler

            # TODO - get HMC acceptance rate here?
            # Is the current method weird, also?

            rate = mh_sampler_acceptance_rate(mh_sampler)
            println("$iteration, $(round(rate, digits=2))")

            mh_sampler.acceptions = 0
            mh_sampler.rejections = 0
        else
            println("$iteration")
        end
    end
end