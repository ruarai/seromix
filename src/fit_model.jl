include("dependencies.jl")

# data_code = ARGS[1]
data_code = "sim_study_hanam_2018_2"

run_dir = "runs/$(data_code)/"

model_data = load("runs/$data_code/model_data.hdf5")
obs_df = DataFrame(model_data["observations"])

p = read_model_parameters(model_data)

maximum(length.(make_obs_views(obs_df)))

model = waning_model(
    p,

    make_obs_lookup(obs_df),
    make_obs_views(obs_df),
    obs_df.observed_titre
);

gibbs_sampler = make_gibbs_sampler(model, :infections, 0.004, p.n_t_steps, p.n_subjects)

chain = @time sample(
    model, gibbs_sampler, 
    MCMCThreads(), 200, 6,

    callback = log_callback
);

plot(chain, [:mu_long, :mu_sum], seriestype = :traceplot)
# plot(chain[1500:end], [:mu_long, :mu_sum], seriestype = :traceplot)

plot(chain, [:sigma_long, :sigma_short], seriestype = :traceplot)
plot(chain, [:tau], seriestype = :traceplot)


save_draws(chain, "$run_dir/chain.parquet")


# This PPD should go to its own function ?
# ppd_obs = make_ppd(chain[1900:end], 1, p)


# fit_data = Dict(
#     "chain" => df_to_tuple(DataFrame(chain)),
#     "ppd" => df_to_tuple(ppd_obs)
# );

# save("$run_dir/fit_data.hdf5", fit_data)