include("dependencies.jl")

rng = Random.Xoshiro(1)

run_dir = "runs/hanam_2018/"
model_data = load("$run_dir/model_data.hdf5")

obs_df = DataFrame(model_data["observations"])

p = read_model_parameters(model_data)

birth_data = DataFrame(model_data["subject_birth_data"])

p.subject_birth_ix .= convert.(Int, birth_data.year_of_birth) .- 1967

prior_infection_dist = MatrixBetaBernoulli(1.0, 1.0, p)
proposal_function = proposal_original_corrected

initial_params = make_initial_params_kucharski_data_study(p, 8, model_data["initial_infections_manual"], rng)
turing_model = waning_model_kucharski

# initial_params = make_initial_params_age(p, obs_df, 4, rng)

model = make_waning_model(p, obs_df; prior_infection_dist = prior_infection_dist, turing_model = turing_model);

gibbs_sampler = make_gibbs_sampler(model, p, proposal_function);

chain = sample_chain(
    model, initial_params, gibbs_sampler, p, rng;
    n_sample = 200_000, n_thinning = 100, n_chain = 8
);

using Plots


ix_start = 1
plot(chain[ix_start:end], [:mu_long, :mu_short], seriestype = :traceplot)
plot(chain_sum_infections(chain, p))

plot(chain[ix_start:end], [:sigma_long, :sigma_short], seriestype = :traceplot)
plot(chain[ix_start:end], [:obs_sd], seriestype = :traceplot)
plot(chain[ix_start:end], [:omega], seriestype = :traceplot)
plot(chain[ix_start:end], [:beta], seriestype = :traceplot)

plot(chain[ix_start:end], [:lp], seriestype = :traceplot)


heatmap(chain_infections_prob(chain[1800:2000], p)')


@gif for i in 1:20:2000
    heatmap(chain_infections_prob(chain[i], p)')
end

chain_name = "loo_test"
save_draws(chain, "$run_dir/chain_$chain_name.parquet")

