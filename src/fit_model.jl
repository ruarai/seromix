include("dependencies.jl")

include("experimental/infer_distances.jl")

data_code = "hanam_2018"
rng = Random.Xoshiro(1)

run_dir = "runs/$(data_code)/"
model_data = load("$run_dir/model_data.hdf5")

obs_df = DataFrame(model_data["observations"])

p = read_model_parameters(model_data)

prior_infection_dist = MatrixBetaBernoulli(1.3, 8.0, p.n_t_steps, p.n_subjects)
proposal_function = propose_swaps_original_corrected!
initial_params = make_initial_params_infer_dist(p, obs_df, 6, rng)

turing_model = waning_model_free

model = make_waning_model(p, obs_df; prior_infection_dist = prior_infection_dist, turing_model = turing_model);

# model = model_full | (sigma_long = 0.13, sigma_short = 0.03);

param_symbols = model_symbols_apart_from(model, [:infections])

gibbs_sampler = Gibbs(
    :infections => make_mh_infection_sampler(p.n_t_steps, p.n_subjects, proposal_function),
    param_symbols => ESS()
);
    

chain = sample_chain(
    model, initial_params, gibbs_sampler, rng;
    n_sample = 500, n_thinning = 1, n_chain = 6
);

ix_start = 1
heatmap(chain_infections_prob(chain[ix_start:end], p)')


plot(chain[ix_start:end], [:mu_long, :mu_short], seriestype = :traceplot)
plot(chain[ix_start:end], [:sigma_long, :sigma_short], seriestype = :traceplot)
plot(chain[ix_start:end], [:obs_sd], seriestype = :traceplot)
plot(chain[ix_start:end], [:omega], seriestype = :traceplot)
plot(chain[ix_start:end], [:tau], seriestype = :traceplot)
plot(chain[ix_start:end], [:dist_scale], seriestype = :traceplot)



scatter(chain[ix_start:end,Symbol("mu_long"),:], chain[ix_start:end,Symbol("strain_locations[2]"),:])

plot(chain, [Symbol("strain_gaps[$i]") for i in 1:5], seriestype = :traceplot)

plot(chain_sum_infections(chain, p))


@gif for i in 1:20:2000
    heatmap(chain_infections_prob(chain[i], p)')
end

chain_name = "infer_distances_1"
save_draws(chain, "$run_dir/chain_$chain_name.parquet")