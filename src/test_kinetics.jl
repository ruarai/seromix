include("dependencies.jl")


t_steps = 1:1:50


n_strain = 3
n_t_steps = length(t_steps)


mu_long = 3.0
mu_short = 5.0
omega = 0.05
sigma_long = 0.2
sigma_short = 0.2
tau = 0.3

dist_matrix = [abs(i - j) for i in 1:n_strain, j in 1:n_strain]

infections = zeros(Bool, n_t_steps, n_strain)
infections[3, 1] = true
infections[15, 3] = true
infections[30, 2] = true

inf_vector = findall(infections)

y = waning_curve(
    mu_long,
    mu_short, omega,
    sigma_long, sigma_short,
    tau,


    dist_matrix,
    inf_vector, n_strain, 
    t_steps
)

plot(y')

plot(y[:, 10])
plot!(y[:, 15])
plot!(y[:, 20])
plot!(y[:, 25])

obs_df = DataFrame(t = Float64[], ix_strain = Int[], titre = Float64[])

for ix_t in 1:3:50, ix_strain in 1:n_strain
    titre = y[ix_strain, ix_t] + rand(Normal(0, 0.3))
    push!(obs_df, (t = t_steps[ix_t], ix_strain = ix_strain, titre = titre))
end

obs_df.ix_t = coalesce.(indexin(obs_df.t, unique(obs_df.t)), -1)

write_parquet("data/obs.parquet", obs_df)

model = waning_model(
    dist_matrix,
    n_strain, n_t_steps,
    
    obs_df, obs_df.titre
)

symbols = DynamicPPL.syms(DynamicPPL.VarInfo(model))
symbols = symbols[findall(symbols .!= :infections)]


mh_sampler = externalsampler(MHInfectionSampler(), adtype=Turing.DEFAULT_ADTYPE, unconstrained=false)

sampler = Gibbs(
    :infections => RepeatSampler(mh_sampler, 2),
    symbols => HMC(0.01, 10)
)


chain = sample(model, sampler, MCMCThreads(), 6000, 4)
save_draws(chain[2000:end], "data/chain.parquet")

plot(chain[:mu_long])
plot(chain[:mu_short])

# plot(chain)

ppd_obs = DataFrame(ix_t = Int[], t = Float64[], ix_strain = Int[])

for ix_t in eachindex(t_steps), ix_strain in 1:n_strain
    push!(ppd_obs, (ix_t = ix_t, t = t_steps[ix_t], ix_strain = ix_strain))
end

# Posterior predictive
pred_model = waning_model(
    dist_matrix,
    n_strain, n_t_steps,
    
    ppd_obs, missing
)

chain_thinning = sample(2000:length(chain), 1000)
ppd = StatsBase.predict(pred_model, chain[chain_thinning])
save_draws(ppd, "data/ppd.parquet")
write_parquet("data/ppd_obs.parquet", ppd_obs)


