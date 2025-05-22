
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


@model function test_titre_model_2(y, t_min, t_max)
    mu ~ Uniform(0, 8)
    sd ~ LogNormal(0, 1)

    for i in eachindex(y)

        y[i] ~ TitreNormal(
            mu, sd,
            convert(typeof(mu), t_min),
            convert(typeof(mu), t_max)
        )
    end

end

dist = TitreArrayNormal(fill(3.75, 10000), 0.5, 0.0, 8.0)
y = rand(dist)
model = test_titre_model_2(y, 0.0, 8.0)

chain = sample(model, NUTS(1000, 0.6), 500);

plot(chain, [:mu, :sd], seriestype = :traceplot)



beta_dist = Beta(1.3, 8)
mean(beta_dist)

plot(0:0.01:1, pdf.(beta_dist, 0:0.01:1))

dist = MatrixBetaBernoulli(1.3, 8.0, 10, 10)

y = rand(dist)
heatmap(y)

ys = [rand(dist) for i in 1:10000]


y_lpdf = [logpdf(dist, y) for y in ys]

y_sums = [sum(y) for y in ys]

histogram(y_sums, bins = 0:1:100)
heatmap(y)