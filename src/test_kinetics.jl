include("dependencies.jl")


t_steps = 0:1:30


n_strain = 5


mu_long = 3.0
mu_short = 5.0
omega = 0.05
sigma_long = 0.2
sigma_short = 0.2
tau = 0.3

dist_matrix = [abs(i - j) for i in 1:n_strain, j in 1:n_strain]

infections = [(1.0, 1), (10.0, 3), (15.0, 5)]

y = waning_curve(
    mu_long,
    mu_short, omega,
    sigma_long, sigma_short,
    tau,


    dist_matrix,
    infections, n_strain, 
    t_steps
)

plot(y[:, 10])
plot!(y[:, 15])
plot!(y[:, 20])
plot!(y[:, 25])

obs_df = DataFrame(t = Float64[], ix_strain = Int[], titre = Float64[])

for ix_t in 1:1:30
    for ix_strain in 1:n_strain
        titre = y[ix_strain, ix_t] + rand(Normal(0, 0.5))
        push!(obs_df, (t = t_steps[ix_t], ix_strain = ix_strain, titre = titre))
    end
end

obs_df.ix_t = coalesce.(indexin(obs_df.t, unique(obs_df.t)), -1)

write_parquet("data/obs.parquet", obs_df)

model = waning_model(
    dist_matrix,
    infections,
    n_strain,
    
    obs_df, obs_df.titre
)

s = rand(model)
logjoint(model, s)


logjoint(
    model,
    (mu_long = 3.0, mu_short = 5.0, omega = 0.05, sigma_long = 0.2, sigma_short = 0.2, tau = 0.3)
)

chain = sample(model, NUTS(1000, 0.65), 1000)
save_draws(chain, "data/chain.parquet")

plot(chain)

ppd_obs = DataFrame(ix_t = Int[], t = Float64[], ix_strain = Int[])

for ix_t in eachindex(t_steps), ix_strain in 1:n_strain
    push!(ppd_obs, (ix_t = ix_t, t = t_steps[ix_t], ix_strain = ix_strain))
end

# Posterior predictive
pred_model = waning_model(
    dist_matrix,
    infections,
    n_strain,
    
    ppd_obs, missing
)

chain_thinning = sample(1:length(chain), 100)
ppd = StatsBase.predict(pred_model, chain[chain_thinning])
save_draws(ppd, "data/ppd.parquet")
write_parquet("data/ppd_obs.parquet", ppd_obs)


