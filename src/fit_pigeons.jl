include("dependencies.jl")

using Pigeons, Bijectors
using Plots

# include("likelihood_model/waning_model_tempering.jl")
include("pigeons/explorer.jl")

data_code = "hanam_2018"
rng = Random.Xoshiro(1)

run_dir = "runs/$(data_code)/"
model_data = load("$run_dir/model_data.hdf5")

obs_df = DataFrame(model_data["observations"])

p = read_model_parameters(model_data)

prior_infection_dist = MatrixBetaBernoulli(1.0, 1.0, p.n_t_steps, p.n_subjects)

model = make_waning_model(
    p, obs_df; prior_infection_dist = prior_infection_dist,
    use_corrected_titre = false,
    # turing_model = waning_model_tempering
);

pigeon_model = TuringLogPotential(model);
const PigeonModelType = typeof(pigeon_model);

function Pigeons.initialization(target::PigeonModelType, rng::AbstractRNG, ::Int64)
    result = DynamicPPL.VarInfo(rng, target.model, DynamicPPL.SampleFromPrior(), DynamicPPL.PriorContext())
    result = DynamicPPL.link(result, target.model)

    Pigeons.update_state!(result, :mu_long, 1, 2.0)
    Pigeons.update_state!(result, :mu_short, 1, 2.0)
    Pigeons.update_state!(result, :omega, 1, 0.8)
    Pigeons.update_state!(result, :sigma_long, 1, 0.1)
    Pigeons.update_state!(result, :sigma_short, 1, 0.05)
    Pigeons.update_state!(result, :tau, 1, 0.05)
    Pigeons.update_state!(result, :obs_sd, 1, 1.5)

    return result
end

symbols_not_inf = model_symbols_apart_from(model, [:infections])

explorer = StateExplorer(
    SliceSampler(),
    proposal_original_corrected,
    [i for i in symbols_not_inf], p.n_t_steps, p.n_subjects, 0.1, 1.0
)

pt = pigeons(
    target = TuringLogPotential(model),
    # n_rounds = 9, n_chains = 16, multithreaded = true,
    n_rounds = 14, n_chains = 64, multithreaded = true,   
    explorer = explorer,
    record = [traces, round_trip, Pigeons.timing_extrema, Pigeons.allocation_extrema]
);

plot(pt.shared.tempering.communication_barriers.localbarrier)

chain = Chains(pt);
plot(chain, [:mu_long, :mu_short], seriestype = :traceplot)
plot(chain, Symbol("infections[20]"), seriestype = :traceplot)

plot(chain, [:sigma_long, :sigma_short], seriestype = :traceplot)
plot(chain, [:obs_sd], seriestype = :traceplot)
plot(chain, [:omega], seriestype = :traceplot)
plot(chain, [:tau], seriestype = :traceplot)


plot(chain, [:log_density], seriestype = :traceplot)


heatmap(chain_infections_prob_2(chain, p)')
# heatmap(model_data["infections_matrix"]')

n_inf = chain_sum_infections(chain, p)

scatter(chain[:mu_long], n_inf)
scatter(chain[:mu_long], chain[:tau])
scatter(chain[:tau], n_inf)
scatter(chain[:mu_short], chain[:omega])

using StatsPlots
@df pt.shared.reports.swap_prs StatsPlots.plot(:round, :mean, group = :first)


chain_name = "pigeons_4"
save_draws(chain, "$run_dir/chain_$chain_name.parquet")
JLD2.save("$run_dir/pt_$chain_name.jld2", Dict("pt" => pt))