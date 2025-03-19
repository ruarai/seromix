

include("dependencies.jl")

data_code = "hanam_2018"

runs_dir = "runs/sim_inf_$(data_code)/"
mkpath(runs_dir)

model_data = load("data/$data_code/model_data.hdf5")

p = read_model_parameters(model_data)

n_t_steps = p.n_t_steps
n_subjects = p.n_subjects

mu_long = 2.0
mu_short = 2.7
omega = 0.8
sigma_long = 0.13
sigma_short = 0.03
tau = 0.04


infections = rand(Bernoulli(0.2), (n_t_steps, n_subjects))

for ix_subject in 1:n_subjects
    if p.subject_birth_ix[ix_subject] > 0
        infections[1:p.subject_birth_ix[ix_subject], ix_subject] .= false
    end
end

heatmap(infections')


infections_df = DataFrame(stack([[i[1], i[2]] for i in findall(infections)])', :auto)
write_parquet("$runs_dir/infections_df.parquet", infections_df)

complete_obs = expand_grid(
    ix_t_obs = 1:n_t_steps, ix_strain = 1:n_t_steps, ix_subject = 1:n_subjects,
    observed_titre = 0.0
)

# obs_lookup::Vector{Dict{Int, Vector{Tuple{Int,Int}}}},
# obs_views::Vector{UnitRange{Int}},
# y::AbstractArray{T}


waning_curve!(
    mu_long, mu_short, omega,
    sigma_long, sigma_short, tau,

    p.antigenic_distances, p.time_diff_matrix, p.subject_birth_ix,

    infections,

    make_obs_lookup(complete_obs), make_obs_views(complete_obs),
    complete_obs.observed_titre
)

write_parquet("$runs_dir/complete_obs.parquet", complete_obs)

function filt_age(ix_subject, ix_t_obs)
    return ix_t_obs > p.subject_birth_ix[ix_subject]
end
    

modelled_years = model_data["modelled_years"]
obs_df = filter(:ix_t_obs => ix_t_obs -> modelled_years[ix_t_obs] % 4 == 0, complete_obs)
obs_df = filter([:ix_subject, :ix_t_obs] => filt_age, obs_df)
#obs_df.observed_titre = round.(max.(0, obs_df.observed_titre .+ rand(Normal(0, 1.3), nrow(obs_df))))
obs_df.observed_titre = obs_df.observed_titre .+ rand(Normal(0, 1.3), nrow(obs_df))

write_parquet("$runs_dir/obs.parquet", obs_df)

model = waning_model(
    p,

    make_obs_lookup(obs_df),
    make_obs_views(obs_df), nrow(obs_df),
    obs_df.observed_titre
);

symbols_not_inf = model_symbols_apart_from(model, :infections)

# HMC step size must be tuned to balance efficiency of the two samplers
gibbs_sampler = Gibbs(
    :infections => make_mh_infection_sampler(),
    symbols_not_inf => HMC(0.006, 10) # Must be reduced with number of individuals?
)

chain = @time sample(
    model, gibbs_sampler, 
    MCMCThreads(), 100, 6,
    callback = log_callback
);

# Save as JLD2?

plot(chain, [:mu_long, :mu_sum], seriestype = :traceplot)
plot(chain, [:sigma_long, :sigma_short], seriestype = :traceplot)
plot(chain, [:tau], seriestype = :traceplot)

save_draws(chain, "$runs_dir/chain.parquet")


ppd_obs = make_ppd(chain[4000:end], 50, p)

write_parquet("$runs_dir/ppd_obs.parquet", ppd_obs)
