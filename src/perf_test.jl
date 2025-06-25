include("dependencies.jl")

using Pigeons, Bijectors
using Plots

# include("likelihood_model/waning_model_tempering.jl")
include("pigeons/explorer.jl")

data_code = "hanam_2018"
rng = Random.Xoshiro(1)

run_dir = "runs/$(data_code)/"
model_data = load("$run_dir/model_data.hdf5")

obs_df = DataFrame(model_data["observations"])

p = read_model_parameters(model_data)

prior_infection_dist = MatrixBetaBernoulli(1.0, 1.0, p.n_t_steps, p.n_subjects)

model = make_waning_model(p, obs_df; prior_infection_dist = prior_infection_dist);

using BenchmarkTools
using LoopVectorization

y_latent = rand(Normal(1.0, 1.5), 100)
x = TitreArrayNormal(y_latent, 1.5, const_titre_min, const_titre_max)
y = sort(rand(x))


apply_logpdf(y, y_latent, 1.5, const_titre_min, const_titre_max)
apply_logpdf_simd(y, y_latent, 1.5, const_titre_min, const_titre_max)


ys = [rand(x) for i in 1:10000]
plot(
    [apply_logpdf(y, y_latent, 1.5, const_titre_min, const_titre_max) for y in ys] .-
    [apply_logpdf_simd(y, y_latent, 1.5, const_titre_min, const_titre_max) for y in ys]
)


@benchmark apply_logpdf(y, $y_latent, 1.5, const_titre_min, const_titre_max) setup=(y = rand(x))
@benchmark apply_logpdf_simd(y, $y_latent, 1.5, const_titre_min, const_titre_max) setup=(y = rand(x))

s = rand(model)
logjoint(model, s)

