

include("dependencies.jl")

t_steps = 1:10
n_t_steps = length(t_steps)
n_ind = 100

mu_long = 2.0
mu_short = 2.7
omega = 0.8
sigma_long = 0.13
sigma_short = 0.03

tau = 0.04

dist_matrix = [abs(i - j) for i in 1:1.0:n_t_steps, j in 1:1.0:n_t_steps]

infections = rand(Bernoulli(0.2), (n_t_steps, n_ind))

infections_df = DataFrame(stack([[i[1], i[2]] for i in findall(infections)])', :auto)
write_parquet("data/infections_df.parquet", infections_df)

complete_obs = expand_grid(t = 1:n_t_steps, s = 1:n_t_steps, i = 1:n_ind)
complete_obs.y = zeros(nrow(complete_obs))

waning_curve!(
    mu_long, mu_short, omega,
    sigma_long, sigma_short, tau,

    dist_matrix,
    infections,

    make_obs_lookup_individuals(complete_obs),
    make_obs_views(complete_obs),
    complete_obs.y
)

write_parquet("data/complete_obs.parquet", complete_obs)


obs_df = filter(:t => t -> t % 1 == 0, complete_obs)
obs_df.y = round.(max.(0, obs_df.y .+ rand(Normal(0, 0.5), nrow(obs_df))))
write_parquet("data/obs.parquet", obs_df)

model = waning_model(
    dist_matrix,

    n_ind, n_t_steps,
    
    make_obs_lookup_individuals(obs_df),
    make_obs_views(obs_df), nrow(obs_df),
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

chain = @time sample(
    model, gibbs_sampler, 
    MCMCThreads(), 2000, 6,
    callback = log_callback
);

# chain = sample(
#     model, gibbs_sampler, 
#     MCMCThreads(), 5000, 8,
#     #num_warmup = 10000,
#     thinning = 10,
#     callback = log_callback
# );

plot(chain, [:mu_long, :mu_sum], seriestype = :traceplot)
plot(chain, [:sigma_long, :sigma_short], seriestype = :traceplot)
plot(chain, [:tau, :omega], seriestype = :traceplot)

save_draws(chain, "data/chain.parquet")


ppd_start_draw_ix = 1000
# Only first 10 individuals for PPD.
ppd_obs = expand_grid(t = 1:n_t_steps, s = 1:n_t_steps, i = 1:10)
# Posterior predictive
pred_model = waning_model(
    dist_matrix,
    n_ind, n_t_steps,
    
    make_obs_lookup_individuals(obs_df),
    make_obs_views(obs_df), nrow(obs_df),
    missing
)

chain_thinning = sample(ppd_start_draw_ix:length(chain), 100);
ppd = StatsBase.predict(pred_model, chain[chain_thinning]);
save_draws(ppd, "data/ppd.parquet")
write_parquet("data/ppd_obs.parquet", ppd_obs)


