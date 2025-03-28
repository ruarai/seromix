include("dependencies.jl")

data_code = ARGS[1]
# data_code = "sim_study_simple_1"

run_dir = "runs/$(data_code)/"

model_data = load("runs/$data_code/model_data.hdf5")
obs_df = DataFrame(model_data["observations"])

p = read_model_parameters(model_data)

n_max_ind_obs = maximum(length.(make_obs_views(obs_df)))

model = waning_model(
    p,

    make_obs_lookup(obs_df),
    make_obs_views(obs_df),
    n_max_ind_obs,
    obs_df.observed_titre
);

gibbs_sampler = make_gibbs_sampler(model, :infections, 0.004, p.n_t_steps, p.n_subjects)

n_thinning = 1
n_sample = 300
chain = @time sample(
    model, gibbs_sampler, 
    MCMCThreads(), n_sample ÷ n_thinning, 6,

    thinning = n_thinning,
    callback = log_callback
);

# plot(chain, [:mu_long, :mu_sum], seriestype = :traceplot)
# plot(chain[1500:end], [:mu_long, :mu_sum], seriestype = :traceplot)

# plot(chain, [:sigma_long, :sigma_short], seriestype = :traceplot)
# plot(chain, [:tau], seriestype = :traceplot)


save_draws(chain, "$run_dir/chain.parquet")