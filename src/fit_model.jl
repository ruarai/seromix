include("dependencies.jl")


data_code = "hanam_2018"
rng = Random.Xoshiro(1)

run_dir = "runs/$(data_code)/"
model_data = load("$run_dir/model_data.hdf5")

obs_df = DataFrame(model_data["observations"])

p = read_model_parameters(model_data)

prior_infection_dist = MatrixBetaBernoulli(1.3, 8.0, p.n_t_steps, p.n_subjects)
proposal_function = propose_swaps_original_corrected!
initial_params = make_initial_params_sim_study(p, obs_df, 6, rng)

turing_model = waning_model_linear

model = make_waning_model(p, obs_df; prior_infection_dist = prior_infection_dist, turing_model = turing_model);

gibbs_sampler = make_gibbs_sampler(model, p, proposal_function)

chain = sample_chain(
    model, initial_params, gibbs_sampler, rng;
    n_sample = 10_000, n_thinning = 5, n_chain = 6
);

# heatmap(model_data["infections_matrix"]')
heatmap(chain_infections_prob(chain[1800:end], p)')

@gif for i in 1:20:2000
    heatmap(chain_infections_prob(chain[i], p)')
end

plot(chain, [:mu_long, :mu_short], seriestype = :traceplot)
plot(chain, [:sigma_long, :sigma_short], seriestype = :traceplot)
plot(chain, [:obs_sd], seriestype = :traceplot)
plot(chain, [:omega], seriestype = :traceplot)
plot(chain, [:tau], seriestype = :traceplot)

plot(chain_sum_infections(chain, p))

chain_name = "linear_basic_2"
save_draws(chain, "$run_dir/chain_$chain_name.parquet")