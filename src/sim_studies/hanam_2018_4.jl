

real_model_data = load("runs/hanam_2018_age/model_data.hdf5")
sp = read_static_parameters(real_model_data)

sim_scenarios = expand_grid(
    ix_sim = 1:3,
    endemic_mean_ar = x_endemic_mean_ar
)
sim_scenarios.ix_row = 1:nrow(sim_scenarios)

for sim_row in eachrow(sim_scenarios)
    run_name = "sim_study_hanam_2018_4/$(sim_row.ix_row)"

    run_dir = mkpath("runs/$(run_name)/")

    attack_rates = vcat(
        rand(rng, LogitNormal(logit(pandemic_mean_ar), sd_ar)),
        [
            rand(rng, LogitNormal(logit(sim_row.endemic_mean_ar), sd_ar))
            for i in 2:sp.n_t_steps
        ]
    )

    infections = infections_from_attack_rate(attack_rates,sp.n_subjects)
    mask_infections_birth_year!(infections,sp.subject_birth_ix) 

    complete_obs = simulate_latent_titre(continuous_params, sp, infections)

    # Only include observations which are available in the study data
    real_obs = DataFrame(real_model_data["observations"])
    observations = innerjoin(complete_obs, real_obs[!, [:ix_t_obs, :ix_strain, :ix_subject]], on = [:ix_t_obs, :ix_strain, :ix_subject])

    apply_titre_obs_noise(observations, rng, continuous_params.obs_sd)

    simulation_metadata = (
        pandemic_mean_ar = pandemic_mean_ar,
        endemic_mean_ar = sim_row.endemic_mean_ar
    )

    model_data = collect_model_data(
        complete_obs, observations,
        infections, continuous_params, 
        real_model_data, simulation_metadata
    )

    save("$run_dir/model_data.hdf5", model_data)
end