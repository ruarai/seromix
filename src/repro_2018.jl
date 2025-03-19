include("dependencies.jl")

data_code = "hanam_2018"

runs_dir = "runs/repro_2018_$(data_code)/"
mkpath(runs_dir)

model_data = load("data/$data_code/model_data.hdf5")
obs_df = DataFrame(model_data["observations_df"])

p = read_model_parameters(model_data)

model = waning_model(
    p,

    make_obs_lookup(obs_df),
    make_obs_views(obs_df), nrow(obs_df),
    obs_df.observed_titre
);

gibbs_sampler = make_gibbs_sampler(model, :infections, 0.008)

chain = @time sample(model, gibbs_sampler, 500);

chain = @time sample(
    model, gibbs_sampler, 
    MCMCThreads(), 5000, 6,
    callback = log_callback
);

plot(chain, [:mu_long, :mu_sum], seriestype = :traceplot)
plot(chain, [:sigma_long, :sigma_short], seriestype = :traceplot)
plot(chain, [:tau], seriestype = :traceplot)

save_draws(chain, "$runs_dir/chain.parquet")

ppd_obs = make_ppd(chain[3000:end], 50, p)

write_parquet("$runs_dir/ppd_obs.parquet", ppd_obs)
