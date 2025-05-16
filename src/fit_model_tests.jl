include("dependencies.jl")

using Plots


# data_code = "sim_study_simple_2"
data_code = "sim_study_hanam_2018_3"

run_dir = "runs/$(data_code)/"

model_data = load("$run_dir/model_data.hdf5")

obs_df = DataFrame(model_data["observations"])

p = read_model_parameters(model_data)

model = make_waning_model(p, obs_df);

# model = model | (infections = Matrix{Bool}(model_data["infections_matrix"]), );
# model = model | (mu_long = 2.0, mu_short = 2.0, sigma_long = 0.15, sigma_short = 0.05, tau = 0.05);

symbols_not_inf = model_symbols_apart_from(model, :infections)


gibbs_sampler = Gibbs(
    :infections => make_mh_infection_sampler(p.n_t_steps, p.n_subjects),
    symbols_not_inf => make_mh_parameter_sampler()
)

# gibbs_sampler = make_gibbs_sampler(model, :infections, 0.005, 10, p)

# infections_0 = initial_infections_matrix(p, obs_df, Random.default_rng())
# chain = @time sample(model, gibbs_sampler, 2000, thinning = 1, initial_params = vec(infections_0));

chain = @time sample_chain(
    model, gibbs_sampler;
    n_sample = 1000, n_thinning = 1, n_chain = 6
);


heatmap(model_data["infections_matrix"]')
# heatmap(infections_0')
heatmap(chain_infections_prob(chain[1000], p)')
heatmap(chain_infections_prob(chain[800:end], p)')

@gif for i in 1:50:5000
    heatmap(chain_infections_prob(chain[i], p)')
end


plot(chain, [:mu_long, :mu_short], seriestype = :traceplot)
plot(chain, [:sigma_long, :sigma_short], seriestype = :traceplot)
plot(chain, [:tau], seriestype = :traceplot)


ix_strain_obs = unique(obs_df.ix_strain)
X = chain_sum_infections(chain, p, 1:500:5000, ix_strain_obs)
plot(X)
hline!([sum(model_data["infections_matrix"][ix_strain_obs, :])])
