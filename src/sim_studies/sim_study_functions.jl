

function collect_infections_df(infections)
    infections_df = DataFrame(stack([[i[1], i[2]] for i in findall(infections)])', :auto)
    rename!(infections_df, ["ix_t", "ix_subject"] )
    return infections_df
end

function simulate_latent_titre(continuous_params, model_parameters, infections)
    # Must be ordered by ix_subject
    complete_obs = expand_grid(
        ix_t_obs = 1:model_parameters.n_t_steps,
        ix_strain = 1:model_parameters.n_t_steps,
        ix_subject = 1:model_parameters.n_subjects,
        observed_titre = 0.0
    )

    waning_curve!(
        continuous_params.mu_long, continuous_params.mu_short, continuous_params.omega,
        continuous_params.sigma_long, continuous_params.sigma_short, continuous_params.tau,

        model_parameters.antigenic_distances,
        model_parameters.time_diff_matrix,
        model_parameters.subject_birth_ix,

        infections,

        # TODO fix with new lookup scheme
        make_obs_lookup(complete_obs), make_obs_views(complete_obs),
        complete_obs.observed_titre
    )

    return complete_obs
end

function infections_from_attack_rate(attack_rates, n_subjects)
    return Matrix(stack([rand(rng, Bernoulli(a), (n_subjects)) for a in attack_rates])')
end

function apply_titre_obs_noise(observations, rng, obs_sd)
    observations.observed_titre = rand(
        rng,
        TitreArrayNormal(
            observations.observed_titre,
            obs_sd, 
            const_titre_min, const_titre_max
        )
    )
end

function collect_model_data(
    complete_obs,
    observations,
    infections,
    continuous_params,
    real_model_data,
    simulation_metadata
)

    infections_df = collect_infections_df(infections)

    return Dict(
        "modelled_years" => real_model_data["modelled_years"],
        "antigenic_distances" => real_model_data["antigenic_distances"],
        "observations" => df_to_tuple(observations),
        "complete_obs" => df_to_tuple(complete_obs),
        "infections" => df_to_tuple(infections_df),
        "infections_matrix" => Matrix{Float64}(infections),
        "subject_birth_data" => real_model_data["subject_birth_data"],
        "simulation_metadata" => [simulation_metadata],
        "continuous_params" => [continuous_params]
    )
end

_round(x) = round(x; digits = 2)

describe_prior_dist(dist::MatrixBernoulli) = "bernoulli_$(_round(dist.p))"
describe_prior_dist(dist::MatrixBetaBernoulli) = "beta_bernoulli_$(_round(dist.alpha))_$(_round(dist.beta))"