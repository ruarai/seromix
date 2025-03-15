include("dependencies.jl")

t_steps = 1:16
n_t_steps = length(t_steps)
n_ind = 100

mu_long = 3.0
mu_short = 5.0
omega = 0.05
sigma_long = 0.2
sigma_short = 0.2

tau = 0.3

dist_matrix = [abs(i - j) for i in 1:1.0:n_t_steps, j in 1:1.0:n_t_steps]

infections = zeros(Bool, n_t_steps, n_ind)

infections = rand(Bernoulli(0.2), (n_t_steps, n_ind))

infections_df = DataFrame(stack([[i[1], i[2]] for i in findall(infections)])', :auto)
write_parquet("data/infections_df.parquet", infections_df)

complete_obs = expand_grid(t = 1:n_t_steps, s = 1:n_t_steps, i = 1:n_ind)

matrix_obs = Matrix(complete_obs[:,1:3])
complete_obs.y .= waning_curve_optimised(
    mu_long, mu_short, omega,
    sigma_long, sigma_short, tau,

    dist_matrix,
    infections,

    make_obs_matrix(complete_obs),
    make_obs_lookup(complete_obs), nrow(complete_obs)
)

write_parquet("data/complete_obs.parquet", complete_obs)


obs_df = filter(:t => t -> t % 4 == 0, complete_obs)
obs_df.y = obs_df.y .+ rand(Normal(0, 0.3), nrow(obs_df))
write_parquet("data/obs.parquet", obs_df)

model = waning_model(
    dist_matrix,

    n_ind, n_t_steps,
    
    make_obs_matrix(obs_df),
    make_obs_lookup(obs_df), nrow(obs_df),
    obs_df.y
)

# using TuringBenchmarking
# benchmark_model(model; check=false, adbackends=[:forwarddiff])

symbols = DynamicPPL.syms(DynamicPPL.VarInfo(model))
symbols = symbols[findall(symbols .!= :infections)]


mh_sampler = externalsampler(MHInfectionSampler(), adtype=Turing.DEFAULT_ADTYPE, unconstrained=false)

# Must somehow balance the level of exploration of the MH sampler
# with that of the HMC sampler -- so repeating MH or changing HMC
# step size
gibbs_sampler = Gibbs(
    :infections => mh_sampler,
    symbols => HMC(0.005, 10) # Must be reduced with number of individuals?
)


# @profview sample(model, gibbs_sampler, 6000);
chain = sample(model, gibbs_sampler, 1000)

chain = sample(model, gibbs_sampler, MCMCThreads(), 12000, 6)

plot(chain, [:mu_long, :mu_sum], seriestype = :traceplot)
plot(chain, [:sigma_long, :sigma_short], seriestype = :traceplot)
plot(chain, [ :tau], seriestype = :traceplot)


skip_n_draws = 6000
save_draws(chain[skip_n_draws:end], "data/chain.parquet")

plot(chain[skip_n_draws:end], [:mu_long, :mu_sum], seriestype = :traceplot)

ppd_obs = expand_grid(t = 1:n_t_steps, s = 1:n_t_steps, i = 1:n_ind)
# Posterior predictive
pred_model = waning_model(
    dist_matrix,
    n_ind, n_t_steps,
    
    make_obs_matrix(ppd_obs),
    make_obs_lookup(ppd_obs), nrow(ppd_obs),
    missing
)

chain_thinning = sample(skip_n_draws:length(chain), 100)
ppd = StatsBase.predict(pred_model, chain[chain_thinning])
save_draws(ppd, "data/ppd.parquet")
write_parquet("data/ppd_obs.parquet", ppd_obs)


