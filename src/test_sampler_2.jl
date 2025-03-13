include("dependencies.jl")


@model function test_model(y_obs)
    p ~ Uniform(0, 1)
    y ~ filldist(Bernoulli(p), length(y_obs))

    for i in eachindex(y_obs)
        y_obs[i] ~ y[i] ? Normal(1, 0.3) : Normal(0, 0.3)
    end
end

y_obs = rand(Bernoulli(0.2), 5) .+ rand(Normal(0, 0.3), 5)
model = test_model(y_obs)

mh_sampler = externalsampler(MHInfectionSampler(), adtype=Turing.DEFAULT_ADTYPE, unconstrained=false)

sampler = Gibbs(
    :y => mh_sampler,
    :p => HMC(0.05, 10)
)

chain = sample(model, sampler, 1000)

plot(chain)
