

function simulate_hanam_2018(
    endemic_mean_ar;
    drop_age = false,
    seed = 1
)
    rng = Random.Xoshiro(seed)

    pandemic_mean_ar = 0.5
    sd_ar = 0.5

    continuous_params = make_continuous_params()

    # As template for amount of data available + age data
    real_model_data = load("runs/hanam_2018_age/model_data.hdf5")
   sp = read_static_parameters(real_model_data)

    attack_rates = vcat(
        rand(rng, LogitNormal(logit(pandemic_mean_ar), sd_ar)),
        [
            rand(rng, LogitNormal(logit(endemic_mean_ar), sd_ar))
            for i in 2:sp.n_t_steps
        ]
    )

    infections = infections_from_attack_rate(rng, attack_rates,sp.n_subjects)
    mask_infections_birth_year!(infections,sp.subject_birth_ix) 
    
    complete_obs = simulate_latent_titre(continuous_params, sp, infections)

    real_obs = DataFrame(real_model_data["observations"])
    observations = innerjoin(complete_obs, real_obs[!, [:ix_t_obs, :ix_strain, :ix_subject]], on = [:ix_t_obs, :ix_strain, :ix_subject])

    apply_titre_obs_noise!(observations, rng, continuous_params.obs_sd)

    model_data = collect_model_data(
        complete_obs, observations,
        infections, continuous_params, 
        real_model_data
    )

    if drop_age
        model_data["subject_birth_data"].ix_t_birth .= 0
    end

    model_data["attack_rates"] = DataFrame(ix_t = 1:sp.n_t_steps, attack_rate = attack_rates)

    return model_data
end

