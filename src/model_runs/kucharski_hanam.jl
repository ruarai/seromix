include("../dependencies.jl")

# Reproduces the HaNam data study from Kucharski (2018)

data_code = "hanam_2018"
rng = Random.Xoshiro(1)

run_dir = "runs/$(data_code)/"
model_data = load("$run_dir/model_data.hdf5")

obs_df = DataFrame(model_data["observations"])

p = read_model_parameters(model_data)

# Kucharski (2018) used an implicit Bernoulli(0.5) prior over infections
# prior_infection_dist = MatrixBernoulli(0.5, p.n_t_steps, p.n_subjects)
prior_infection_dist = MatrixBetaBernoulli(1.3, 8.0, p.n_t_steps, p.n_subjects)
# and a proposal function which omitted the hastings ratio
# proposal_function = propose_swaps_original_no_hastings_ratio!
proposal_function = propose_swaps_original_corrected!
# and initial infections from sim study.
initial_params = make_initial_params_data_study(6, model_data["initial_infections_manual"], rng)

model = make_waning_model(p, obs_df; prior_infection_dist = prior_infection_dist);
gibbs_sampler = make_gibbs_sampler(model, p, proposal_function)

# Sample 50,000 draws with a thinning factor of 25, with 6 chains each
chain = sample_chain(
    model, initial_params, gibbs_sampler, rng;
    n_sample = 50000, n_thinning = 25, n_chain = 6
);

chain_name = "prior_beta_1.3_8.0_corrected"

save_draws(chain, "$run_dir/chain_$chain_name.parquet")