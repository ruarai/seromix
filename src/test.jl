
using Distributions
using Turing

using Plots, StatsPlots


y_obs = rand(Normal(1.5, 0.5), 100)

@model function test_model(y_obs)
    mu ~ Normal(0, 1.0)

    for i in eachindex(y_obs)
        
        y_obs[i] ~ Normal(mu, 0.5)
    end
end

model = test_model(y_obs)

sampler = HMC(0.05, 10)

chain = sample(model, sampler, 1000)




plot(chain)
