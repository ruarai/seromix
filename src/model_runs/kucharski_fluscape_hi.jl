include("../dependencies.jl")

# Reproduces the fluscape (HI) data study from Kucharski (2018)

data_code = "fluscape_2009_HI"
rng = Random.Xoshiro(1)

run_dir = "runs/$(data_code)/"
model_data = load("$run_dir/model_data.hdf5")

obs_df = DataFrame(model_data["observations"])

p = read_model_parameters(model_data)

# Kucharski (2018) used an implicit Bernoulli(0.5) prior over infections
prior_infection_dist = MatrixBernoulli(0.5, p.n_t_steps, p.n_subjects)
# and a proposal function which omitted the hastings ratio
proposal_function = proposal_original_uncorrected
# and initial values with some variance (but infections initialised as in sim study)
initial_params = make_initial_params_kucharski_data_study_fluscape(p, obs_df, 6, rng)

model_free = make_waning_model(p, obs_df; prior_infection_dist = prior_infection_dist);

# Fix certain parameters (note in original study, there is between-chain variance in omega)
model = model_free | (
    mu_short = 1e-10, omega = 0.8
) 

gibbs_sampler = make_gibbs_sampler(model, p, proposal_function)

chain = sample_chain(
    model, initial_params, gibbs_sampler, rng;
    n_sample = 50000, n_thinning = 25, n_chain = 6
);


chain_name = "prior_50_uncorrected"

save_draws(chain, "$run_dir/chain_$chain_name.parquet")