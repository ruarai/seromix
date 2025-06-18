include("dependencies.jl")

using Pigeons, Bijectors
using Plots

include("likelihood_model/waning_model_tempering.jl")

data_code = "hanam_2018"
rng = Random.Xoshiro(1)

run_dir = "runs/$(data_code)/"
model_data = load("$run_dir/model_data.hdf5")

obs_df = DataFrame(model_data["observations"])

p = read_model_parameters(model_data)

prior_infection_dist = MatrixBetaBernoulli(1.0, 1.0, p.n_t_steps, p.n_subjects)

model = make_waning_model(
    p, obs_df; prior_infection_dist = prior_infection_dist,
    # turing_model = waning_model_tempering
);

pigeon_model = TuringLogPotential(model);
const PigeonModelType = typeof(pigeon_model);

# TODO remove?
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

pt = pigeons(
    target = TuringLogPotential(model),
    n_rounds = 12, n_chains = 32, multithreaded = true,
    # n_rounds = 1, n_chains = 1, multithreaded = false,
    explorer = CustomExplorer2(SliceSampler(1.0, 20, 3, 4096), symbols_not_inf, p.n_t_steps, p.n_subjects),
    record = [traces]
);

pt = increment_n_rounds!(pt, 1)
pt = pigeons(pt)


plot(pt.shared.tempering.communication_barriers.localbarrier)


chain = Chains(pt);
ix_start = 1
plot(chain[ix_start:end], [:mu_long, :mu_short], seriestype = :traceplot)
plot(chain[ix_start:end], Symbol("infections[20]"), seriestype = :traceplot)

plot(chain[ix_start:end], [:sigma_long, :sigma_short], seriestype = :traceplot)
plot(chain[ix_start:end], [:obs_sd], seriestype = :traceplot)
plot(chain[ix_start:end], [:omega], seriestype = :traceplot)
plot(chain[ix_start:end], [:tau], seriestype = :traceplot)


plot(chain[ix_start:end], [:log_density], seriestype = :traceplot)


inf_mat = chain_infections_prob_2(chain, p)
heatmap(inf_mat')

n_inf = chain_sum_infections(chain, p)
plot(n_inf)


scatter(chain[:mu_long], n_inf)
scatter(chain[:mu_long], chain[:tau])
scatter(chain[:tau], n_inf)
scatter(chain[:mu_short], chain[:omega])

using StatsPlots
@df pt.shared.reports.swap_prs StatsPlots.plot(:round, :mean, group = :first, legend = false)


chain_name = "pigeons_1"
save_draws(chain, "$run_dir/chain_$chain_name.parquet")