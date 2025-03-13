
include("dependencies.jl")

@model function test_model(y_obs)
    y ~ filldist(Bernoulli(0.2), length(y_obs))

    sigma ~ LogNormal(0, 1)

    for i in eachindex(y_obs)
        if y[i]
            y_obs[i] ~ Normal(1, sigma)
        else
            y_obs[i] ~ Normal(0, sigma)
        end
    end
end

y_true = rand(Bernoulli(0.2), 10)

y_obs = y_true .+ rand(Normal(0, 0.1), length(y_true))
model = test_model(y_obs)


function prop_fn(rng, x)
    return rand(rng, Bernoulli(0.5), length(x))
end

f_prop = FunctionProposal{true}(prop_fn)

gibbs = Gibbs(
    :y => RepeatSampler(MH(:y => f_prop), 10),
    :sigma => NUTS(2000, 0.65)
)

chain = sample(model, gibbs, MCMCThreads(), 8000, 4)
plot(chain[6000:end][:sigma])
plot(chain[6000:end])