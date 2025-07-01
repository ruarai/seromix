include("dependencies.jl")


data_code = "hanam_2018"
rng = Random.Xoshiro(1)

run_dir = "runs/$(data_code)/"
model_data = load("$run_dir/model_data.hdf5")

obs_df = DataFrame(model_data["observations"])

p = read_model_parameters(model_data)

prior_infection_dist = MatrixBetaBernoulliTimeVarying(1.0, 1.0, p)
proposal_function = proposal_original_corrected

initial_params = make_initial_params_kucharski_data_study(p, 4, model_data["initial_infections_manual"], rng)

model = make_waning_model(p, obs_df; prior_infection_dist = prior_infection_dist);

gibbs_sampler = make_gibbs_sampler(model, p, proposal_function);

symbols_not_inf = model_symbols_apart_from(model, [:infections])
gibbs_sampler = Gibbs(
    :infections => make_mh_infection_sampler(p, proposal_function),
    symbols_not_inf => externalsampler(RobustAdaptiveMetropolis2(0.234, 0.6, nothing, 0.0, Inf))
)

chain = sample_chain(
    model, initial_params, gibbs_sampler, p, rng;
    n_sample = 100000, n_thinning = 50, n_chain = 4
);

heatmap(chain_infections_prob(chain[1800:2000], p)')

using Plots

ix_start = 1
plot(chain[ix_start:end], [:mu_long, :mu_short], seriestype = :traceplot)
plot(chain[ix_start:end], [:lp], seriestype = :traceplot)
plot(chain[ix_start:end], [:sigma_long, :sigma_short], seriestype = :traceplot)
plot(chain[ix_start:end], [:obs_sd], seriestype = :traceplot)
plot(chain[ix_start:end], [:omega], seriestype = :traceplot)
plot(chain[ix_start:end], [:tau], seriestype = :traceplot)


plot(chain_sum_infections(chain, p))


@gif for i in 1:20:2000
    heatmap(chain_infections_prob(chain[i], p)')
end

# chain_name = "prior_beta_1.3_8.0_corrected"
# save_draws(chain, "$run_dir/chain_$chain_name.parquet")