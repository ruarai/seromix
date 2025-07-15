include("dependencies.jl")

rng = Random.Xoshiro(1)

run_dir = "runs/hanam_2018/"
model_data = load("$run_dir/model_data.hdf5")

obs_df = DataFrame(model_data["observations"])

sp = read_static_parameters(model_data)

birth_data = DataFrame(model_data["subject_birth_data"])

prior_infection_dist = MatrixBetaBernoulli(1.0, 1.0, sp)
proposal_function = proposal_original_corrected

# prior_infection_dist = MatrixBernoulli(0.5, sp)
# proposal_function = proposal_original_uncorrected


turing_model = waning_model_kucharski
# initial_params = make_initial_params_kucharski_data_study(sp, 4, rng, model_data["initial_infections_manual"])
# initial_params = make_initial_params_broad(sp, 4, rng)
initial_params = make_initial_params_kucharski_sim_study(sp, 4, rng, obs_df)

# turing_model = waning_model_age_effect
# initial_params = make_initial_params_age(sp, 16, rng, obs_df)

model = make_waning_model(
   sp, obs_df; prior_infection_dist = prior_infection_dist, turing_model = turing_model,
   mixture_importance_sampling = false
);

gibbs_sampler = make_gibbs_sampler(model, sp, proposal_function);
# gibbs_sampler = make_gibbs_sampler_original(model, sp, proposal_function);
# gibbs_sampler = make_gibbs_sampler_slice(model, sp, proposal_function);

# symbols_not_inf = model_symbols_apart_from(model, [:infections])

# gibbs_sampler = Gibbs(
#     :infections => make_mh_infection_sampler(sp, proposal_function; prop_sample = 1.0, n_repeats = 5),
#     symbols_not_inf => NUTS(100, 0.65)
# )

chain = sample_chain(
    model, initial_params, gibbs_sampler, sp, rng;
    n_sample = 20000, n_thinning = 10, n_chain = 4
);


set_lp!(model, chain)

chain = chain[:,:,15]

using Plots
ix_start = 1
plot(chain[ix_start:end], [:mu_long, :mu_short], seriestype = :traceplot)
plot(chain_sum_infections(chain, sp))

plot(chain[ix_start:end], [:sigma_long, :sigma_short], seriestype = :traceplot)
plot(chain[ix_start:end], [:obs_sd], seriestype = :traceplot)
plot(chain[ix_start:end], [:omega], seriestype = :traceplot)
plot(chain[ix_start:end], [:tau], seriestype = :traceplot)

plot(chain[ix_start:end], [:lp], seriestype = :traceplot)

heatmap(chain_infections_prob(chain[1800:2000], sp)')


@gif for i in 1:20:2000
    heatmap(chain_infections_prob(chain[i], sp)')
end

# chain_name = "hmc_test"
# save_draws(chain, "$run_dir/chain_$chain_name.parquet")




# chain_df = DataFrame(chain[1000:end])
# lpp = model_sum_mixIS(chain_df, sp, obs_df, model)



