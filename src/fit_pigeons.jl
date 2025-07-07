include("dependencies.jl")

using Pigeons, Bijectors
using Plots

include("pigeons/explorer.jl")

data_code = "hanam_2018_age"
rng = Random.Xoshiro(1)

run_dir = "runs/$(data_code)/"
model_data = load("$run_dir/model_data.hdf5")

obs_df = DataFrame(model_data["observations"])

p = read_model_parameters(model_data)

prior_infection_dist = MatrixBetaBernoulli(1.0, 1.0, p)

turing_model = waning_model_kucharski

model = make_waning_model(
    p, obs_df; prior_infection_dist = prior_infection_dist, turing_model = turing_model
);

pt_target = TuringLogPotential(model)

function Pigeons.initialization(target::typeof(pt_target), rng::AbstractRNG, ::Int64)
    result = DynamicPPL.VarInfo(rng, target.model, DynamicPPL.SampleFromPrior(), DynamicPPL.PriorContext())

    Pigeons.update_state!(result, :mu_long, 1, 2.0)
    Pigeons.update_state!(result, :mu_short, 1, 2.0)
    Pigeons.update_state!(result, :omega, 1, 0.8)
    Pigeons.update_state!(result, :sigma_long, 1, 0.1)
    Pigeons.update_state!(result, :sigma_short, 1, 0.05)
    Pigeons.update_state!(result, :tau, 1, 0.05)
    Pigeons.update_state!(result, :obs_sd, 1, 1.5)

    return DynamicPPL.link(result, target.model)
end

# Approx timing
# sum([0.5 * 2^i for i in 1:15]) / (60 * 60)

symbols_not_inf = model_symbols_apart_from(model, [:infections])
pt = pigeons(
    target = pt_target,
    n_rounds = 12, n_chains = 16, multithreaded = true,   
    explorer = GibbsExplorer(proposal_original_corrected, [i for i in symbols_not_inf], p),
    extended_traces = true,
    record = [traces, round_trip, Pigeons.timing_extrema, Pigeons.allocation_extrema, index_process]
);


chain = Chains(pt);


# chain_name = "pigeons_5"
# save_draws(chain, "$run_dir/chain_$chain_name.parquet")
# JLD2.save("$run_dir/pt_$chain_name.jld2", Dict("pt" => pt))