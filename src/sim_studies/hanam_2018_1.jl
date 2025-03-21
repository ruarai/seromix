
include("../dependencies.jl")

data_code = "hanam_2018_1"

run_dir = "runs/sim_study_$(data_code)/"
mkpath(run_dir)

model_data = load("data/hanam_2018/model_data.hdf5")

p = read_model_parameters(model_data)

n_t_steps = p.n_t_steps
n_subjects = p.n_subjects

mu_long = 2.0
mu_short = 2.0
omega = 0.75
sigma_long = 0.2
sigma_short = 0.1
tau = 0.05

infections = rand(Bernoulli(0.2), (n_t_steps, n_subjects))

for ix_subject in 1:n_subjects
    if p.subject_birth_ix[ix_subject] > 0
        infections[1:p.subject_birth_ix[ix_subject], ix_subject] .= false
    end
end

heatmap(infections')


infections_df = DataFrame(stack([[i[1], i[2]] for i in findall(infections)])', :auto)
rename!(infections_df, ["t", "ix_subject"] )
write_parquet("$run_dir/infections_df.parquet", infections_df)

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

write_parquet("$run_dir/complete_obs.parquet", complete_obs)

function filt_age(ix_subject, ix_t_obs)
    return ix_t_obs > p.subject_birth_ix[ix_subject]
end
    

modelled_years = model_data["modelled_years"]
obs_df = filter(:ix_t_obs => ix_t_obs -> modelled_years[ix_t_obs] % 4 == 0, complete_obs)
obs_df = filter([:ix_subject, :ix_t_obs] => filt_age, obs_df)
#obs_df.observed_titre = round.(max.(0, obs_df.observed_titre .+ rand(Normal(0, 1.3), nrow(obs_df))))
obs_df.observed_titre = obs_df.observed_titre .+ rand(Normal(0, 1.3), nrow(obs_df))

write_parquet("$run_dir/obs.parquet", obs_df)