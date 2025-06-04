include("dependencies.jl")


data_code = "sim_study_hanam_2018_1"
rng = Random.Xoshiro(1)

run_dir = "runs/$(data_code)/"
model_data = load("$run_dir/model_data.hdf5")

obs_df = DataFrame(model_data["observations"])

p = read_model_parameters(model_data)


prior_infection_dist = MatrixBetaBernoulli(1.3, 8.0, p.n_t_steps, p.n_subjects)
# prior_infection_dist = MatrixBernoulli(0.15, p.n_t_steps, p.n_subjects)

proposal_function = proposal_jitter

initial_params = make_initial_params_kucharski_sim_study(p, obs_df, 8, rng)

model = make_waning_model(p, obs_df; prior_infection_dist = prior_infection_dist);


symbols_not_inf = model_symbols_apart_from(model, [:infections])

gibbs_sampler = Gibbs(
    :infections => make_mh_infection_sampler(p.n_t_steps, p.n_subjects, proposal_function),
    symbols_not_inf => HMC(0.005, 10)
)


chain = sample_chain(
    model, initial_params, gibbs_sampler, rng;
    n_sample = 2000, n_thinning = 1, n_chain = 8
);

heatmap(chain_infections_prob(chain[1800:end], p)')


ix_start = 1
plot(chain[ix_start:end], [:mu_long, :mu_short], seriestype = :traceplot)
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