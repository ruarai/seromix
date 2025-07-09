include("dependencies.jl")

rng = Random.Xoshiro(1)

run_dir = "runs/sim_study_tiny_1/"
model_data = load("$run_dir/model_data.hdf5")

obs_df = DataFrame(model_data["observations"])

sp = read_fixed_parameters(model_data)

birth_data = DataFrame(model_data["subject_birth_data"])

prior_infection_dist = MatrixBetaBernoulli(1.0, 1.0, sp)
proposal_function = proposal_original_corrected

turing_model = waning_model_age_effect
initial_params = make_initial_params_age(sp, obs_df, 8, rng)

# turing_model = waning_model_kucharski
# initial_params = make_initial_params_kucharski_data_study(sp, 4, model_data["initial_infections_manual"], rng)

model = make_waning_model(
   sp, obs_df; prior_infection_dist = prior_infection_dist, turing_model = turing_model,
   mixture_importance_sampling = false
);

gibbs_sampler = make_gibbs_sampler(model, sp, proposal_function);

chain = sample_chain(
    model, initial_params, gibbs_sampler, sp, rng;
    n_sample = 20_000, n_thinning = 10, n_chain = 8
    # n_sample = 100, n_thinning = 1, n_chain = 8
);

set_lp!(model, chain)

using Plots
ix_start = 1
plot(chain[ix_start:end], [:mu_long, :mu_short], seriestype = :traceplot)
plot(chain_sum_infections(chain, sp))

plot(chain[ix_start:end], [:sigma_long, :sigma_short], seriestype = :traceplot)
plot(chain[ix_start:end], [:obs_sd], seriestype = :traceplot)
plot(chain[ix_start:end], [:omega], seriestype = :traceplot)
plot(chain[ix_start:end], [:tau, :beta, :intercept], seriestype = :traceplot)

plot(chain[ix_start:end], [:lp], seriestype = :traceplot)


heatmap(chain_infections_prob(chain[1800:2000], sp)')


@gif for i in 1:20:2000
    heatmap(chain_infections_prob(chain[i], sp)')
end

# chain_name = "nonlinear_test"
# save_draws(chain, "$run_dir/chain_$chain_name.parquet")




chain_df = DataFrame(chain[1500:end])
lpp = model_sum_mixIS(chain_df, sp, obs_df, model)