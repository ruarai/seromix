

include("dependencies.jl")


n_t_steps = 30
n_strain = 5

infections = [(1, 1), (10, 3), (15, 5)]

y_data = zeros(n_t_steps, n_strain)
infection_binary = zeros(n_t_steps, n_strain)

for ix_inf in eachindex(infections)
    t_inf, strain = infections[ix_inf] 
    infection_binary[t_inf, strain] = true

    for t in t_inf:n_t_steps
        y_data[t, strain] = 1.0
    end
end
y_data = y_data .+ rand(Normal(0.0, 0.1), size(y_data))


heatmap(infection_binary')
heatmap(y_data')


@model function infection_model(
    y_obs,

    n_t_steps, n_strain
)
    p_inf ~ Uniform(0, 1)

    infection_binary ~ filldist(Bernoulli(p_inf), n_t_steps, n_strain)

    for ix_strain in 1:n_strain
        was_inf = false

        for ix_t in 1:n_t_steps
            was_inf = was_inf || infection_binary[ix_t, ix_strain]

            if was_inf
                y_obs[ix_t, ix_strain] ~ Normal(1.0, 0.1)
            else
                y_obs[ix_t, ix_strain] ~ Normal(0.0, 0.1)
            end
        end
    end
end


model = infection_model(
    y_data, 

    n_t_steps, n_strain
)

x = rand(model)
logjoint(model, x)


gibbs_sampler = Gibbs(HMC(0.2, 3, :p_inf), PG(20, :infection_binary))

sampler = PG(1000)

# Need to create a custom sampler here which 
# can do the special proposals over the infection history
# Then have NUTS/HMC do the rest.

# sampler = MH()

chain = sample(model, gibbs_sampler, 200);
save_draws(chain, "data/chain.parquet")

pred_model = infection_model(
    fill(missing, size(y_data)), 

    n_t_steps, n_strain
)

chain_thinning = sample(100:length(chain), 100)
ppd = StatsBase.predict(pred_model, chain[chain_thinning])

save_draws(ppd, "data/ppd.parquet")
