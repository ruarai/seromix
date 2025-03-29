include("../dependencies.jl")

data_code = "hanam_2018_3"

run_dir = "runs/sim_study_$(data_code)/"
mkpath(run_dir)

real_model_data = load("runs/hanam_2018/model_data.hdf5")

p = read_model_parameters(real_model_data)

n_t_steps = p.n_t_steps
n_subjects = p.n_subjects

mu_long = 2.0
mu_short = 2.0
omega = 0.75
sigma_long = 0.2
sigma_short = 0.1
tau = 0.05
sigma_obs = 1.5


modelled_years = real_model_data["modelled_years"]

# sd_param = 0.5
# mean_offset = (sd_param / 2) ^ 2 / 2

# attack_rates = vcat(
#     rand(LogNormal(log(0.5) - mean_offset, sd_param)),
#     rand(LogNormal(log(0.15) - mean_offset, sd_param), length(modelled_years) - 1)
# )
# attack_rates = fill(0.2, length(modelled_years))
# infections = Matrix(stack([rand(Bernoulli(a), (n_subjects)) for a in attack_rates])')
# for ix_subject in 1:n_subjects
#     if p.subject_birth_ix[ix_subject] > 0
#         infections[1:p.subject_birth_ix[ix_subject], ix_subject] .= false
#     end
# end

infections = rand(Bernoulli(0.2), (n_t_steps, n_subjects))
mask_infections_birth_year!(infections, p.subject_birth_ix)
heatmap(infections')


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

    make_obs_lookup(complete_obs), make_obs_views(complete_obs),
    complete_obs.observed_titre
)

function filt_age(ix_subject, ix_t_obs)
    return ix_t_obs > p.subject_birth_ix[ix_subject]
end
    

observed_strains = unique(DataFrame(real_model_data["observations"]).ix_strain)


real_observations = DataFrame(real_model_data["observations"])

# Only include observations that were in the real study.
observations = innerjoin(
    complete_obs, 
    real_observations[:, [:ix_subject, :ix_strain, :ix_t_obs]],
    on = [:ix_subject, :ix_strain, :ix_t_obs]
)

observations.observed_titre = rand(TitreArrayNormal(observations.observed_titre, sigma_obs, const_titre_min, const_titre_max))

model_data = Dict(
    "modelled_years" => modelled_years,
    "antigenic_distances" => real_model_data["antigenic_distances"],
    "observations" => df_to_tuple(observations),
    "complete_obs" => df_to_tuple(complete_obs),
    "infections" => df_to_tuple(infections_df),
    "subject_birth_data" => real_model_data["subject_birth_data"]
)

save("$run_dir/model_data.hdf5", model_data)