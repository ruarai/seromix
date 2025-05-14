include("dependencies.jl")

using Plots


data_code = "sim_study_simple_2"

run_dir = "runs/$(data_code)/"

model_data = load("$run_dir/model_data.hdf5")

obs_df = DataFrame(model_data["observations"])

p = read_model_parameters(model_data)

model = make_waning_model(p, obs_df);

symbols_not_inf = model_symbols_apart_from(model, :infections)

gibbs_sampler = Gibbs(
    :infections => make_mh_infection_sampler(p.n_t_steps, p.n_subjects),
    symbols_not_inf => make_mh_parameter_sampler()
)

# chain = @time sample(model, gibbs_sampler, 4000);

chain = @time sample_chain(
    model, gibbs_sampler;
    n_sample = 1000, n_thinning = 4, n_chain = 6
);

heatmap(model_data["infections_matrix"]')
heatmap(chain_infections_prob(chain[200:end], p)')

@gif for i in 1:5:1000
    heatmap(chain_infections_prob(chain[i], p)')
end


plot(chain, [:mu_long, :mu_sum], seriestype = :traceplot)
plot(chain, [:sigma_long, :sigma_short], seriestype = :traceplot)
plot(chain, [:tau], seriestype = :traceplot)

plot(chain, [:mu_long, :mu_sum], seriestype = :traceplot)

plot(chain_sum_infections(chain)) # TODO this is broken?
hline!([sum(model_data["infections_matrix"])])