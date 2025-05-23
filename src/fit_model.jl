include("dependencies.jl")


data_code = "fluscape_2009_neuts"
rng = Random.Xoshiro(1)

run_dir = "runs/$(data_code)/"
model_data = load("$run_dir/model_data.hdf5")

obs_df = DataFrame(model_data["observations"])

p = read_model_parameters(model_data)

prior_infection_dist = MatrixBernoulli(0.5, p.n_t_steps, p.n_subjects)
proposal_function = propose_swaps_original_no_hastings_ratio!
initial_params = make_initial_params_sim_study(p, obs_df, 6, rng)

model = make_waning_model(p, obs_df; prior_infection_dist = prior_infection_dist);
gibbs_sampler = make_gibbs_sampler(model, p, proposal_function)

chain = sample_chain(
    model, initial_params, gibbs_sampler, rng;
    n_sample = 10000, n_thinning = 5, n_chain = 6
);

heatmap(model_data["infections_matrix"]')
heatmap(chain_infections_prob(chain[1800:end], p)')

@gif for i in 1:20:2000
    heatmap(chain_infections_prob(chain[i], p)')
end

plot(chain, [:mu_long, :mu_short], seriestype = :traceplot)
plot(chain, [:sigma_long, :sigma_short], seriestype = :traceplot)
plot(chain, [:tau], seriestype = :traceplot)
plot(chain, [:obs_sd], seriestype = :traceplot)
plot(chain, [:omega], seriestype = :traceplot)


plot(chain_sum_infections(chain, p))
hline!([nrow(model_data["infections"])])

chain_name = "prior_50_uncorrected"

save_draws(chain, "$run_dir/chain_$chain_name.parquet")