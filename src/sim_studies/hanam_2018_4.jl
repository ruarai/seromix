
include("../dependencies.jl")
include("sim_study_functions.jl")

rng = Random.Xoshiro(1)

pandemic_mean_ar = 0.5
x_endemic_mean_ar = 0.1:0.1:0.5
sd_ar = 0.5

continuous_params = (
    mu_long = 2.0,
    mu_short = 2.0,
    omega = 0.75,
    sigma_long = 0.15,
    sigma_short = 0.05,
    tau = 0.05,
    obs_sd = 1.5
)

real_model_data = load("runs/hanam_2018_age/model_data.hdf5")
p = read_model_parameters(real_model_data)

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
            for i in 2:n_t_steps
        ]
    )

    infections = infections_from_attack_rate(attack_rates, p.n_subjects)
    mask_infections_birth_year!(infections, p.subject_birth_ix) 

    complete_obs = simulate_latent_titre(continuous_params, p, infections)

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