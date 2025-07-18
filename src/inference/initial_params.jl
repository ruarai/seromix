

function make_initial_params_kucharski_data_study(sp, n_chain, rng, init_matrix)
    init_matrix_masked = copy(init_matrix)
    mask_infections_birth_year!(init_matrix_masked, sp.subject_birth_ix)

    return [(
        mu_long = 2.0 + rand(rng, Uniform(-0.2, 0.2)),
        mu_short = 2.5 + rand(rng, Uniform(-0.2, 0.2)), 
        omega = 0.8 + rand(rng, Uniform(-0.05, 0.05)), 
        sigma_long = 0.1 + rand(rng, Uniform(-0.02, 0.02)),
        sigma_short = 0.05 + rand(rng, Uniform(-0.005, 0.005)), 
        tau = 0.05 + rand(rng, Uniform(-0.01, 0.01)), 
        obs_sd = 1.5 + rand(rng, Uniform(-0.1, 0.1)), 

        infections = init_matrix_masked
    ) for i in 1:n_chain]
end


function make_initial_params_broad(sp, n_chain, rng)
    return [(
        mu_long = rand(rng, Uniform(0.5, 5.0)),
        mu_short = rand(rng, Uniform(0.5, 5.0)), 
        omega = rand(rng, Uniform(0.5, 1.0)), 
        sigma_long = rand(rng, Uniform(0.0, 0.5)),
        sigma_short = rand(rng, Uniform(0.0, 0.5)), 
        tau = rand(rng, Uniform(0.0, 0.2)), 
        obs_sd = rand(rng, Uniform(3.0, 5.0)), 

        infections = rand(rng, MatrixBetaBernoulli(1.0, 1.0, sp))
    ) for i in 1:n_chain]
end

function make_initial_params_kucharski_sim_study(sp, n_chain, rng, obs_df)
    return [(
        mu_long = 2.0 + rand(rng, Uniform(-0.2, 0.2)),
        mu_short = 2.5 + rand(rng, Uniform(-0.2, 0.2)), 
        omega = 0.8 + rand(rng, Uniform(-0.05, 0.05)), 
        sigma_long = 0.15 + rand(rng, Uniform(-0.02, 0.02)),
        sigma_short = 0.05 + rand(rng, Uniform(-0.005, 0.005)), 
        tau = 0.05 + rand(rng, Uniform(-0.01, 0.01)), 
        obs_sd = 1.5 + rand(rng, Uniform(-0.1, 0.1)), 

        infections = initial_infections_matrix(sp, obs_df, rng)
    ) for i in 1:n_chain]
end


function make_initial_params_kucharski_data_study_fluscape(sp, n_chain, rng, obs_df)
    return [(
        mu_long = 2.0 + rand(rng, Uniform(-0.2, 0.2)),
        tau = 0.05 + rand(rng, Uniform(-0.01, 0.01)), 
        sigma_long = 0.15 + rand(rng, Uniform(-0.02, 0.02)),
        obs_sd = 1.5 + rand(rng, Uniform(-0.1, 0.1)), 

        infections = initial_infections_matrix(sp, obs_df, rng)
    ) for i in 1:n_chain]
end

# Translated from Kucharski model, but not sure that it's necessary.
function initial_infections_matrix(sp, obs_df, rng)
    infections_0 = zeros(Bool,sp.n_t_steps,sp.n_subjects)

    obs_df_grouped = groupby(obs_df, :ix_subject)

    for ix_subject in 1:sp.n_subjects
        
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

    mask_infections_birth_year!(infections_0,sp.subject_birth_ix)

    return infections_0
end


function make_initial_params_age(sp, n_chain, rng, obs_df)
    return [(
        mu_long = 2.0 + rand(rng, Uniform(-0.2, 0.2)),
        mu_short = 2.5 + rand(rng, Uniform(-0.2, 0.2)), 
        omega = 0.8 + rand(rng, Uniform(-0.05, 0.05)), 
        sigma_long = 0.1 + rand(rng, Uniform(-0.02, 0.02)),
        sigma_short = 0.05 + rand(rng, Uniform(-0.005, 0.005)), 
        beta = 0.02 + rand(rng, Uniform(-0.01, 0.01)), 
        tau = 0.05 + rand(rng, Uniform(-0.01, 0.01)), 
        intercept = 0.0 + rand(rng, Uniform(-0.1, 0.1)), 
        obs_sd = 1.5 + rand(rng, Uniform(-0.1, 0.1)), 

        infections = initial_infections_matrix(sp, obs_df, rng)
    ) for i in 1:n_chain]
end