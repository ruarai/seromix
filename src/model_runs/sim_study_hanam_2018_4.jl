include("../dependencies.jl")
include("../sim_studies/sim_study_functions.jl")
rng = Random.Xoshiro(1)


run_dir_base = "sim_study_hanam_2018_4"
name_suffix = "_uncorrected"

run_indices = drop_nothing(tryparse.(Int, readdir("runs/$(run_dir_base)/")))
n_chains = 4

# Assumed same across all scenarios
p = read_model_parameters(load("runs/$(run_dir_base)/1/model_data.hdf5"))

# p.subject_birth_ix .= 0
# println("Setting all birth_ix to zero")

exps = DataFrame(ix_run = Int[], dist = Distribution[])

for ix_run in run_indices
    push!(exps, (ix_run, MatrixBernoulli(0.5, p.n_t_steps, p.n_subjects)))
    push!(exps, (ix_run, MatrixBernoulli(0.3, p.n_t_steps, p.n_subjects)))
    push!(exps, (ix_run, MatrixBernoulli(0.1, p.n_t_steps, p.n_subjects)))

    push!(exps, (ix_run, MatrixBetaBernoulli(1.0, 1.0, p.n_t_steps, p.n_subjects)))
    push!(exps, (ix_run, MatrixBetaBernoulli(1.0, 8.0, p.n_t_steps, p.n_subjects)))
end

n_parallel = Threads.nthreads() ÷ n_chains * 2
println("Expected sampling time of $(round(10 / 60 * ceil(nrow(exps) / n_parallel); digits = 2)) hr")
jobs = split_vector_indices(nrow(exps), n_parallel)

Threads.@threads for row_ix_job in jobs
    println("Running jobs $row_ix_job")
    for exp_row in eachrow(exps[row_ix_job, :])
        run_dir = "runs/$run_dir_base/$(exp_row.ix_run)"
        model_data = load("$run_dir/model_data.hdf5")

        obs_df = DataFrame(model_data["observations"])
        chain_name = describe_prior_dist(exp_row.dist)

        initial_params = make_initial_params_kucharski_sim_study(p, obs_df, n_chains, rng)

        model = make_waning_model(p, obs_df; prior_infection_dist = exp_row.dist);
        gibbs_sampler = make_gibbs_sampler(model, p, proposal_original_uncorrected)

        println("Running $run_dir, $chain_name$name_suffix")

        chain = sample_chain(
            model, initial_params, gibbs_sampler, p, rng;
            n_sample = 100_000, n_thinning = 50, n_chain = n_chains,
            progress = false
        );

        save_draws(chain, "$run_dir/chain_$chain_name$name_suffix.parquet")

        println("Completed $run_dir, $chain_name$name_suffix")
    end
end


