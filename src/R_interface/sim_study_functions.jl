function make_continuous_params()
    return (
        mu_long = 2.0,
        mu_short = 2.0,
        omega = 0.75,
        sigma_long = 0.15,
        sigma_short = 0.05,
        tau = 0.05,
        obs_sd = 1.5
    )
end


function infections_from_attack_rate(rng, attack_rates, n_subjects)
    return Matrix(stack([rand(rng, Bernoulli(a), (n_subjects)) for a in attack_rates])')
end


function simulate_latent_titre(continuous_params, sp, infections; individual_waning_function = individual_waning_kucharski!)
    # Must be ordered by ix_subject
    complete_obs = expand_grid(
        ix_t_obs = 1:sp.n_t_steps,
        ix_strain = 1:sp.n_t_steps,
        ix_subject = 1:sp.n_subjects,
        observed_titre = 0.0
    )

    obs_lookup_strain, obs_lookup_ix = make_obs_lookup(complete_obs)

    waning_curve!(
        continuous_params,
        individual_waning_function,

       sp.antigenic_distances,
       sp.time_diff_matrix,
       sp.subject_birth_ix,

        infections,
        obs_lookup_strain, obs_lookup_ix, make_obs_views(complete_obs),

        complete_obs.observed_titre
    )

    return complete_obs
end

function apply_titre_obs_noise!(observations, rng, obs_sd)
    observations.observed_titre = rand(
        rng,
        TitreArrayNormal(
            observations.observed_titre,
            obs_sd, 
            const_titre_min, const_titre_max
        )
    )
end


function collect_infections_df(infections)
    infections_df = DataFrame(stack([[i[1], i[2]] for i in findall(infections)])', :auto)
    rename!(infections_df, ["ix_t", "ix_subject"] )
    return infections_df
end

function collect_model_data(
    complete_obs,
    observations,
    infections,
    continuous_params,
    real_model_data
)
    infections_df = collect_infections_df(infections)

    return Dict(
        "modelled_years" => real_model_data["modelled_years"],
        "antigenic_distances" => real_model_data["antigenic_distances"],
        "observations" => observations,
        "complete_obs" => complete_obs,
        "infections" => infections_df,
        "infections_matrix" => Matrix{Float64}(infections),
        "subject_birth_data" => DataFrame(real_model_data["subject_birth_data"]),
        "continuous_params" => continuous_params
    )
end