include("dependencies.jl")

# data_code = ARGS[1]
data_code = "hanam_2018"

run_dir = "runs/$(data_code)/"

model_data = load("$run_dir/model_data.hdf5")

obs_df = DataFrame(model_data["observations"])

p = read_model_parameters(model_data)

model = make_waning_model(p, obs_df);

# Condition on true values
# model = model | 
#     (mu_sum = 4.0, mu_long = 2.0, sigma_long = 0.15, sigma_short = 0.05, tau = 0.05);

gibbs_sampler = make_gibbs_sampler(model, :infections, 0.007, 10, p.n_t_steps, p.n_subjects);

chain = @time sample_chain(
    model, gibbs_sampler;
    n_sample = 2000, n_thinning = 1, n_chain = 6
);


heatmap(chain_infections_prob(chain[1500:end], p)')

plot(chain, [:mu_long, :mu_sum], seriestype = :traceplot)
plot(chain, [:sigma_long, :sigma_short], seriestype = :traceplot)
plot(chain, [:tau], seriestype = :traceplot)

save_draws(chain, "$run_dir/chain.parquet")