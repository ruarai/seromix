include("../dependencies.jl")

pandemic_mean_ar = 0.5
x_endemic_ar = 0.1:0.1:0.5
pandemic_sd_ar = 1.0
endemic_sd_ar = 0.5


expand_grid(
    endemic_ar = x_endemic_ar,
    ix_sim = 1:5
)


run_name = "sim_study_hanam_2018_4"
rng = Random.Xoshiro(1)

run_dir = "runs/$(run_name)/"
mkpath(run_dir)

real_model_data = load("runs/hanam_2018/model_data.hdf5")

p = read_model_parameters(real_model_data)

n_t_steps = p.n_t_steps
n_subjects = p.n_subjects

continuous_params = (
    mu_long = 2.0,
    mu_short = 2.0,
    omega = 0.75,
    sigma_long = 0.15,
    sigma_short = 0.05,
    tau = 0.05,
    obs_sd = 1.5
)



modelled_years = real_model_data["modelled_years"]

inf_sd_param = 0.5
mean_offset = (inf_sd_param / 2) ^ 2 / 2

attack_rates = vcat(
    rand(LogitNormal(logit(pandemic_mean_ar), pandemic_sd_ar)),
    [
        rand(LogitNormal(logit(endemic_ar), endemic_sd_ar))
        for i in 2:n_t_steps
    ]
)

infections = Matrix(stack([rand(rng, Bernoulli(a), (n_subjects)) for a in attack_rates])')

mask_infections_birth_year!(infections, p.subject_birth_ix) 
heatmap(infections')


infections_df = DataFrame(stack([[i[1], i[2]] for i in findall(infections)])', :auto)
rename!(infections_df, ["ix_t", "ix_subject"] )

complete_obs = expand_grid(
    ix_t_obs = 1:n_t_steps, ix_strain = 1:n_t_steps, ix_subject = 1:n_subjects,
    observed_titre = 0.0
)

waning_curve!(
    continuous_params.mu_long, continuous_params.mu_short, continuous_params.omega,
    continuous_params.sigma_long, continuous_params.sigma_short, continuous_params.tau,

    p.antigenic_distances, p.time_diff_matrix, p.subject_birth_ix,

    infections,

    make_obs_lookup(complete_obs), make_obs_views(complete_obs),
    complete_obs.observed_titre
)

# Only include observations which are available in the study data
real_obs = DataFrame(real_model_data["observations"])

observations = innerjoin(complete_obs, real_obs[!, [:ix_t_obs, :ix_strain, :ix_subject]], on = [:ix_t_obs, :ix_strain, :ix_subject])

observations.observed_titre = rand(
    rng,
    TitreArrayNormal(observations.observed_titre, obs_sd, const_titre_min, const_titre_max)
)

model_data = Dict(
    "modelled_years" => modelled_years,
    "antigenic_distances" => real_model_data["antigenic_distances"],
    "observations" => df_to_tuple(observations),
    "complete_obs" => df_to_tuple(complete_obs),
    "infections" => df_to_tuple(infections_df),
    "infections_matrix" => Matrix{Float64}(infections),
    "subject_birth_data" => real_model_data["subject_birth_data"]
)

save("$run_dir/model_data.hdf5", model_data)