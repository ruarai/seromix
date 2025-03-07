include("dependencies.jl")



t_steps = 0:1:30

t_infected = [5, 20]

mu = 5.0
omega = 0.03

y = waning_curve(mu, omega, t_infected, t_steps)
y_noisy = y .+ rand(Normal(0.0, 0.5), length(y))

ix_obs = [10, 15, 20, 25]
ix_obs = 1:25
t_obs = t_steps[ix_obs]
titer_obs = y_noisy[ix_obs]

plot(t_steps, y)
scatter!(t_obs, titer_obs)




model = waning_model(t_obs, t_infected, titer_obs)

sampler = NUTS(1000, 0.65)
chain = sample(model, sampler, 1000)

save_draws(chain, "data/chain.parquet")


# Posterior predictive
pred_model = waning_model(t_steps, t_infected, missing)

chain_thin = sample(1:length(chain), 100)
ppd = StatsBase.predict(pred_model, chain[chain_thin])
save_draws(ppd, "data/ppd_df.parquet")


