include("../dependencies.jl")
rng = Random.Xoshiro(1)


run_dir_base = "sim_study_hanam_2018_4"
run_indices = parse.(Int, readdir("runs/$(run_dir_base)/"))

# Assumed same across all scenarios
p = read_model_parameters(load("runs/$(run_dir_base)/1/model_data.hdf5"))

exps = DataFrame(ix_run = Int[], name = String[], dist = Distribution[])

for i in run_indices
    push!(exps, (i, "bernoulli_0.5", MatrixBernoulli(0.5, p.n_t_steps, p.n_subjects)))
    push!(exps, (i, "beta_bernoulli_1.0_1.0", MatrixBetaBernoulli(1.0, 1.0, p.n_t_steps, p.n_subjects)))
    push!(exps, (i, "beta_bernoulli_1.0_8.0", MatrixBetaBernoulli(1.0, 8.0, p.n_t_steps, p.n_subjects)))
end


for exp_row in eachrow(exps)
    run_dir = "runs/$(run_dir_base)/$(exp_row.ix_run)"
    model_data = load("$run_dir/model_data.hdf5")

    obs_df = DataFrame(model_data["observations"])

    prior_infection_dist = exp_row.dist

    proposal_function = propose_swaps_original_corrected!

    initial_params = make_initial_params_sim_study(p, obs_df, 6, rng)

    model = make_waning_model(p, obs_df; prior_infection_dist = prior_infection_dist);
    gibbs_sampler = make_gibbs_sampler(model, p, proposal_function)

    chain = sample_chain(
        model, initial_params, gibbs_sampler, rng;
        n_sample = 20000, n_thinning = 10, n_chain = 6
    );

    save_draws(chain, "$run_dir/chain_$(exp_row.name).parquet")

end