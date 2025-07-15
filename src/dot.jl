




m = 1.5:0.1:3.0

param_symbols = model_symbols_apart_from(model, [:infections, :mu_long, :tau])


chains = 1:4
n_draws = 1980:2000

lp = zeros(length(ms), length(n_draws) * 4)

ix_chain, ix_draw = 4, 2000

draws = DataFrame(chain[:,:,ix_chain])


col_names = names(draws)
ix_infections = findall(s -> startswith(s, "infections"), col_names)
draw = draws[ix_draw,:]

params = NamedTuple{param_symbols}([draw[i] for i in param_symbols])

infections = reshape(Vector(draw[ix_infections]),sp.n_t_steps,sp.n_subjects)
infections = convert(Matrix{Bool}, convert.(Bool, infections))

mu_longs = 1.5:0.01:2.5
taus = 0.02:0.0002:0.06

lp = [logjoint(model, merge(params, (mu_long = mu_long, tau = tau, infections = infections))) for mu_long in mu_longs, tau in taus]



plot(chain[ix_start:end], [:mu_long, :tau], seriestype = :traceplot)

heatmap(taus, mu_longs, lp, color = :viridis)
scatter!(chain[1500:end][:tau], chain[1500:end][:mu_long], mc = "black", ms = 0.5)