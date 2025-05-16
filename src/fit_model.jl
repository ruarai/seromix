include("dependencies.jl")

# data_code = ARGS[1]
data_code = "hanam_2018"

run_dir = "runs/$(data_code)/"

model_data = load("$run_dir/model_data.hdf5")

obs_df = DataFrame(model_data["observations"])

p = read_model_parameters(model_data)

model = make_waning_model(p, obs_df);

gibbs_sampler = make_gibbs_sampler(model, :infections);

chain = @time sample_chain(
    model, gibbs_sampler;
    n_sample = 50000, n_thinning = 25, n_chain = 6
);

heatmap(model_data["infections_matrix"]')
heatmap(chain_infections_prob(chain[1800:end], p)')

@gif for i in 1:20:2000
    heatmap(chain_infections_prob(chain[i], p)')
end

plot(chain, [:mu_long, :mu_short, :sigma_long, :sigma_short, :tau], seriestype = :traceplot)


plot(chain, [:mu_long, :mu_short], seriestype = :traceplot)
plot(chain, [:sigma_long, :sigma_short], seriestype = :traceplot)
plot(chain, [:tau], seriestype = :traceplot)
plot(chain, [:obs_sd], seriestype = :traceplot)
plot(chain, [:omega], seriestype = :traceplot)

plot(chain_sum_infections(chain))
hline!([sum(model_data["infections_matrix"])])

save_draws(chain, "$run_dir/chain.parquet")