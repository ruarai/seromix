
include("dependencies.jl")


@model function test_titre_model(y, t_min, t_max)
    mu ~ Uniform(0, 8)
    sd ~ LogNormal(0, 1)

    y ~ TitreArrayNormal(
        fill(mu, length(y)), sd,
        convert(typeof(mu), t_min),
        convert(typeof(mu), t_max)
    )

end

t_min = 0.0
t_max = 8.0


dist = TitreArrayNormal(fill(7.0, 1), 1.5, t_min, t_max)
sum([exp.(logpdf(dist, [i])) for i in 0:1:8]) # Should be 1



dist = TitreArrayNormal(fill(0.75, 10000), 0.5, t_min, t_max)
y = rand(dist)

model = test_titre_model(y, t_min, t_max)

chain = sample(model, NUTS(1000, 0.6), 500);

plot(chain, [:mu, :sd], seriestype = :traceplot)