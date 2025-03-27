

include("../dependencies.jl")

data_code = "simple_1"

run_dir = "runs/sim_study_$(data_code)/"
mkpath(run_dir)

modelled_years = collect(1980:1999)
n_t_steps = length(modelled_years)
n_subjects = 20

time_diff_matrix = make_time_diff_matrix(modelled_years)

antigenic_distance = abs.(time_diff_matrix)

subject_birth_ix = floor.(Int, reverse(1:n_subjects) .* 0.3)

subject_birth_data = DataFrame(ix_subject = 1:n_subjects, ix_t_birth = subject_birth_ix)

p = FixedModelParameters(
    n_t_steps, n_subjects,
    antigenic_distance,
    time_diff_matrix,
    subject_birth_ix
)

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
rename!(infections_df, ["ix_t", "ix_subject"] )

complete_obs = expand_grid(
    ix_t_obs = 1:n_t_steps, ix_strain = 1:n_t_steps, ix_subject = 1:n_subjects,
    observed_titre = 0.0
)

make_obs_lookup(complete_obs)[1][10]

waning_curve!(
    mu_long, mu_short, omega,
    sigma_long, sigma_short, tau,

    p.antigenic_distances, p.time_diff_matrix, p.subject_birth_ix,

    infections,

    make_obs_lookup(complete_obs), make_obs_views(complete_obs),
    complete_obs.observed_titre
)

function filt_age(ix_subject, ix_t_obs)
    return ix_t_obs >= p.subject_birth_ix[ix_subject]
end
    

observations = filter([:ix_subject, :ix_t_obs] => filt_age, complete_obs)
observations.observed_titre = observations.observed_titre .+ rand(Normal(0, 1.5), nrow(observations))

model_data = Dict(
    "modelled_years" => modelled_years,
    "antigenic_distances" => antigenic_distance,
    "observations" => df_to_tuple(observations),
    "complete_obs" => df_to_tuple(complete_obs),
    "infections" => df_to_tuple(infections_df),
    "subject_birth_data" => df_to_tuple(subject_birth_data)
)

save("$run_dir/model_data.hdf5", model_data)