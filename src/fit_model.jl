include("dependencies.jl")


data_code = "sim_study_simple_hierarchical_1"
rng = Random.Xoshiro(1)

run_dir = "runs/$(data_code)/"
model_data = load("$run_dir/model_data.hdf5")

obs_df = DataFrame(model_data["observations"])

p = read_model_parameters(model_data)

proposal_function = propose_swaps_original_corrected!
initial_params = make_initial_params_sim_study(p, obs_df, 6, rng)

model = make_waning_model(p, obs_df);

symbols_not_inf = model_symbols_apart_from(model, [:infections])

gibbs_sampler = Gibbs(
    :infections => make_mh_infection_sampler(p.n_t_steps, p.n_subjects, proposal_function),
    symbols_not_inf => ESS()
)

chain = sample_chain(
    model, initial_params, gibbs_sampler, rng;
    n_sample = 2000, n_thinning = 1, n_chain = 6
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


plot(chain, [Symbol("row_means[$i]") for i in 1:10], seriestype = :traceplot, ylim = (-4, 2))
plot(chain, [Symbol("row_means[$i]") for i in 20:30], seriestype = :traceplot, ylim = (-4, 2))
plot(model_data["attack_rates_log_odds"][20:30])

plot(chain_sum_infections(chain))
hline!([nrow(model_data["infections"])])

chain_name = "prior_hierarchical"

save_draws(chain, "$run_dir/chain_$chain_name.parquet")