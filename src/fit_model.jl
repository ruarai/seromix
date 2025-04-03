include("dependencies.jl")

data_code = ARGS[1]
# data_code = "sim_study_simple_1"

run_dir = "runs/$(data_code)/"

model_data = load("$run_dir/model_data.hdf5")

obs_df = DataFrame(model_data["observations"])

p = read_model_parameters(model_data)

model = make_waning_model(p, obs_df);

gibbs_sampler = make_gibbs_sampler(model, :infections, 0.0075, p.n_t_steps, p.n_subjects)

chain = @time sample_chain(
    model, gibbs_sampler;
    n_sample = 2500, n_thinning = 10, n_chain = 6
);



plot(chain, [:mu_long, :mu_sum], seriestype = :traceplot)
# plot(chain[1500:end], [:mu_long, :mu_sum], seriestype = :traceplot)

plot(chain, [:sigma_long, :sigma_short], seriestype = :traceplot)
plot(chain, [:tau], seriestype = :traceplot)


save_draws(chain, "$run_dir/chain.parquet")