include("dependencies.jl")

rng = Random.Xoshiro(1)

run_dir = "runs/hanam_2018/"
model_data = load("$run_dir/model_data.hdf5")

obs_df = DataFrame(model_data["observations"])

sp = read_static_parameters(model_data)

birth_data = DataFrame(model_data["subject_birth_data"])

# prior_infection_dist = MatrixBetaBernoulliTimeVarying(1.3, 8.0, sp)

prior_infection_dist = MatrixBetaBernoulli(1.0, 1.0, sp)
proposal_function = proposal_original_corrected

# turing_model = waning_model_kucharski_diff
turing_model = waning_model_kucharski

# turing_model = waning_model_kucharski_priors
# initial_params = make_initial_params_broad(sp, 4, rng)
initial_params = make_initial_params_kucharski_data_study(sp, 4, rng, model_data["initial_infections_manual"])

model = make_waning_model(
   sp, obs_df; prior_infection_dist = prior_infection_dist, turing_model = turing_model,
   mixture_importance_sampling = false
);

gibbs_sampler = make_gibbs_sampler(model, sp, proposal_function);

# symbols_not_inf = model_symbols_apart_from(model, [:infections])

# gibbs_sampler = Gibbs(
    # :infections => make_mh_infection_sampler(sp, proposal_function; prop_sample = 1.0, n_repeats = 5),
    # symbols_not_inf => NUTS(1000, 0.65)
# )

chain = sample_chain(
    model, initial_params, gibbs_sampler, sp, rng;
    n_sample = 2000, n_thinning = 1, n_chain = 4
);



using Plots
ix_start = 1
plot(chain[ix_start:end], [:mu_long, :mu_short], seriestype = :traceplot)
plot(chain_sum_infections(chain, sp))

plot(chain[ix_start:end], [:sigma_long, :sigma_short], seriestype = :traceplot)
plot(chain[ix_start:end], [:obs_sd], seriestype = :traceplot)
plot(chain[ix_start:end], [:omega], seriestype = :traceplot)
plot(chain[ix_start:end], [:tau], seriestype = :traceplot)

plot(chain[ix_start:end], [:lp], seriestype = :traceplot)

heatmap(chain_infections_prob(chain, sp)')

# X = chain_infections_prob(chain[1000:2000], sp)
# # mask_infections_birth_year!(X, sp.subject_birth_ix; mask_val = NaN)
# heatmap(X[:,sortperm(sp.subject_birth_ix)]')


# @gif for i in 1:20:2000
#     heatmap(chain_infections_prob(chain[i], sp)')
# end

chain_sub = chain[1000:end, :, :];
set_lp!(model, chain_sub)

plot(chain_sub, [:lp], seriestype = :traceplot)

scatter(chain_sub[:mu_long], chain_sum_infections(chain_sub, sp))
scatter(chain_sub[:mu_short], chain_sum_infections(chain_sub, sp))

scatter(chain_sub[:mu_long], chain_sum_infections(chain_sub, sp))

plot([sum(chain_infections_prob(chain_sub[:,:,i], sp),dims = 2) for i in 1:16])
plot([sum(chain_infections_prob(chain_sub[:,:,i], sp),dims = 1)' for i in 1:16])


# chain_name = "2025-07-16_1"
# save_draws(chain, "$run_dir/chain_$chain_name.parquet")