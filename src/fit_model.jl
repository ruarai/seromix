include("dependencies.jl")


data_code = "hanam_2018"
rng = Random.Xoshiro(1)

run_dir = "runs/$(data_code)/"
model_data = load("$run_dir/model_data.hdf5")

obs_df = DataFrame(model_data["observations"])

p = read_model_parameters(model_data)


# prior_infection_dist = MatrixBetaBernoulli(1.3, 8.0, p.n_t_steps, p.n_subjects)
prior_infection_dist = MatrixBernoulli(0.15, p.n_t_steps, p.n_subjects)

proposal_function = propose_swaps_original_corrected!

initial_params = make_initial_params_sim_study(p, obs_df, 6, rng)

model = make_waning_model(p, obs_df; prior_infection_dist = prior_infection_dist);
gibbs_sampler = make_gibbs_sampler(model, p, proposal_function)

chain = sample_chain(
    model, initial_params, gibbs_sampler, rng;
    n_sample = 10_000, n_thinning = 5, n_chain = 6
);

ix_start = 1
heatmap(chain_infections_prob(chain[ix_start:end], p)')


plot(chain[ix_start:end], [:mu_long, :mu_short], seriestype = :traceplot)
plot(chain[ix_start:end], [:sigma_long, :sigma_short], seriestype = :traceplot)
plot(chain[ix_start:end], [:obs_sd], seriestype = :traceplot)
plot(chain[ix_start:end], [:omega], seriestype = :traceplot)
plot(chain[ix_start:end], [:tau], seriestype = :traceplot)


plot(chain_sum_infections(chain, p))


@gif for i in 1:20:2000
    heatmap(chain_infections_prob(chain[i], p)')
end

chain_name = "prior_beta_1.3_8.0_corrected"
save_draws(chain, "$run_dir/chain_$chain_name.parquet")