include("dependencies.jl")

using Pigeons, Bijectors
using Plots

include("pigeons/explorer.jl")

data_code = "hanam_2018_age" # NOTE USING AGE here
rng = Random.Xoshiro(1)

run_dir = "runs/$(data_code)/"
model_data = load("$run_dir/model_data.hdf5")

obs_df = DataFrame(model_data["observations"])

sp = read_static_parameters(model_data)

# prior_infection_dist = MatrixBetaBernoulli(1.0, 1.0, sp)
prior_infection_dist = MatrixBetaBernoulliTimeVarying(1.3, 8.0, sp)

turing_model = waning_model_kucharski

model = make_waning_model(
   sp, obs_df; prior_infection_dist = prior_infection_dist, turing_model = turing_model#,
    # mixture_importance_sampling = true
);

pt_target = TuringLogPotential(model)

initial_inf = initial_infections_matrix(sp, obs_df, rng)

function Pigeons.initialization(target::typeof(pt_target), rng::AbstractRNG, ::Int64)
    result = DynamicPPL.VarInfo(rng, target.model, DynamicPPL.SampleFromPrior(), DynamicPPL.PriorContext())

    Pigeons.update_state!(result, :mu_long, 1, 2.0)
    Pigeons.update_state!(result, :mu_short, 1, 2.0)
    Pigeons.update_state!(result, :omega, 1, 0.8)
    Pigeons.update_state!(result, :sigma_long, 1, 0.1)
    Pigeons.update_state!(result, :sigma_short, 1, 0.05)
    Pigeons.update_state!(result, :tau, 1, 0.05)
    Pigeons.update_state!(result, :obs_sd, 1, 1.5)

    for i in 1:(sp.n_t_steps * sp.n_subjects)
        Pigeons.update_state!(result, :infections, i, initial_inf[i])
    end

    return DynamicPPL.link(result, target.model)
end


# Approx timing
sum([2.0 * 2^i for i in 1:11]) / (60 * 60)

symbols_not_inf = model_symbols_apart_from(model, [:infections])
pt = pigeons(
    target = pt_target,
    # n_rounds = 15, n_chains = 64, multithreaded = true,   
    n_rounds = 7, n_chains = 64, multithreaded = true,
    explorer = GibbsExplorer4(proposal_original_corrected, [i for i in symbols_not_inf], 5, sp),
    record = [traces, round_trip, Pigeons.timing_extrema, Pigeons.allocation_extrema, index_process]
);

pt = increment_n_rounds!(pt, 1)
pt = pigeons(pt)


chain_name = "pigeons_5_mixis"
chain = Chains(pt);
save_draws(chain, "$run_dir/chain_$chain_name.parquet")
JLD2.save("$run_dir/pt_$chain_name.jld2", Dict("pt" => pt))