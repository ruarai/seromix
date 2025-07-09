include("../dependencies.jl")

data_code = "sim_study_hanam_2018_3"
rng = Random.Xoshiro(1)

run_dir = "runs/$(data_code)/"
mkpath(run_dir)

real_model_data = load("runs/hanam_2018/model_data.hdf5")

p = read_fixed_parameters(real_model_data)

p.subject_birth_ix .= 0

subject_birth_data = DataFrame(real_model_data["subject_birth_data"])
subject_birth_data.ix_t_birth .= 0
subject_birth_data.year_of_birth .= 1960

n_t_steps = p.n_t_steps
n_subjects = p.n_subjects

mu_long = 2.0
mu_short = 2.0
omega = 0.75
sigma_long = 0.15
sigma_short = 0.05
tau = 0.05
obs_sd = 1.5


modelled_years = real_model_data["modelled_years"]

inf_sd_param = 0.5
mean_offset = (inf_sd_param / 2) ^ 2 / 2

attack_rates = vcat(
    rand(rng, LogNormal(log(0.5) - (inf_sd_param / 2) ^ 2 / 2, inf_sd_param / 2)),
    rand(rng, LogNormal(log(0.15) - inf_sd_param ^ 2 / 2, inf_sd_param), length(modelled_years) - 1)
)
infections = Matrix(stack([rand(rng, Bernoulli(a), (n_subjects)) for a in attack_rates])')

# Disable temporarily to try match Kucharski?
# mask_infections_birth_year!(infections, p.subject_birth_ix) 
# heatmap(infections')

infections_df = DataFrame(stack([[i[1], i[2]] for i in findall(infections)])', :auto)
rename!(infections_df, ["ix_t", "ix_subject"] )

complete_obs = expand_grid(
    ix_t_obs = 1:n_t_steps, ix_strain = 1:n_t_steps, ix_subject = 1:n_subjects,
    observed_titre = 0.0
)

waning_curve!(
    mu_long, mu_short, omega,
    sigma_long, sigma_short, tau,

    p.antigenic_distances, p.time_diff_matrix, p.subject_birth_ix,

    infections,

    # TODO fix with new lookup scheme
    make_obs_lookup(complete_obs), make_obs_views(complete_obs),
    complete_obs.observed_titre
)

observed_strains = unique(DataFrame(real_model_data["observations"]).ix_strain)


# Only include observations from 2007 onwards
observations = filter(:ix_t_obs => ix_t_obs -> modelled_years[ix_t_obs] >= 2007, complete_obs)
observations = filter(:ix_strain => ix_strain -> in(ix_strain, observed_strains), observations)



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
    "subject_birth_data" => df_to_tuple(subject_birth_data)
)

save("$run_dir/model_data.hdf5", model_data)