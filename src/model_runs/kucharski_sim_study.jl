include("../dependencies.jl")

# Reproduces the simulation study from Kucharski (2018)

# Use simulation study data as produced by the original code (to be extra-sure)
data_code = "sim_study_hanam_ak_1"
rng = Random.Xoshiro(1)

run_dir = "runs/$(data_code)/"
model_data = load("$run_dir/model_data.hdf5")

obs_df = DataFrame(model_data["observations"])

p = read_model_parameters(model_data)

# Kucharski (2018) used an implicit Bernoulli(0.5) prior over infections
prior_infection_dist = MatrixBernoulli(0.5, p.n_t_steps, p.n_subjects)
# and a proposal function which omitted the hastings ratio
proposal_function = proposal_original_uncorrected
# and initial values based on observations (I don't think this is too important though)
initial_params = make_initial_params_sim_study(p, obs_df, 6, rng)


model = make_waning_model(p, obs_df; prior_infection_dist = prior_infection_dist);
gibbs_sampler = make_gibbs_sampler(model, p, proposal_function)

# Sample 50,000 draws with a thinning factor of 25, with 6 chains each
chain = sample_chain(
    model, initial_params, gibbs_sampler, rng;
    n_sample = 50000, n_thinning = 25, n_chain = 6
);

chain_name = "prior_50_uncorrected"

save_draws(chain, "$run_dir/chain_$chain_name.parquet")