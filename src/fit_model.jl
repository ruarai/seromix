include("dependencies.jl")

data_code = "sim_study_simple_1"

runs_dir = "runs/$(data_code)/"

model_data = load("runs/$data_code/model_data.hdf5")
obs_df = DataFrame(model_data["observations"])

p = read_model_parameters(model_data)

prior_inf = fill(0.01, (p.n_t_steps, p.n_subjects))

model = waning_model(
    p,

    prior_inf,

    make_obs_lookup(obs_df),
    make_obs_views(obs_df), nrow(obs_df),
    obs_df.observed_titre
);

gibbs_sampler = make_gibbs_sampler(model, :infections, 0.006)

# sample(model, gibbs_sampler, 10, callback = log_callback);
# @profview sample(model, gibbs_sampler, 500, callback = log_callback);

chain = @time sample(
    model, gibbs_sampler, 
    MCMCThreads(), 2000, 6,
    callback = log_callback
);

plot(chain, [:mu_long, :mu_sum], seriestype = :traceplot)
plot(chain, [:sigma_long, :sigma_short], seriestype = :traceplot)
plot(chain, [:tau], seriestype = :traceplot)

save_draws(chain, "$runs_dir/chain.parquet")

ppd_obs = make_ppd(chain[1800:end], 50, p)

write_parquet("$runs_dir/ppd_obs.parquet", ppd_obs)
