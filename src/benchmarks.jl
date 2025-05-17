
include("dependencies.jl")

t_steps = 1:4
n_t_steps = length(t_steps)
n_subjects = 50

mu_long = 2.0
mu_short = 2.0
omega = 0.75
sigma_long = 0.2
sigma_short = 0.1
tau = 0.05

dist_matrix = [abs(i - j) for i in 1:1.0:n_t_steps, j in 1:1.0:n_t_steps]

infections = rand(Bernoulli(0.5), (n_t_steps, n_subjects))

infections_df = DataFrame(stack([[i[1], i[2]] for i in findall(infections)])', :auto)
complete_obs = expand_grid(ix_strain = 1:n_t_steps, ix_t_obs = 1:n_t_steps, ix_subject = 1:n_subjects)

n_obs = nrow(complete_obs)

obs_views = make_obs_views(complete_obs)
obs_lookup = make_obs_lookup(complete_obs)

y = zeros(n_obs)

time_diff_matrix = make_time_diff_matrix(1:n_t_steps)
subject_birth_ix = zeros(Int, n_subjects)

waning_curve!(
    mu_long, mu_short, omega,
    sigma_long, sigma_short, tau,


    dist_matrix, time_diff_matrix,
    subject_birth_ix,

    infections,

    obs_lookup, obs_views,
    y
)

using BenchmarkTools

@benchmark waning_curve_in_place!(
    $mu_long, $mu_short, $omega,
    $sigma_long, $sigma_short, $tau,


    $dist_matrix, $time_diff_matrix,
    $subject_birth_ix,

    $infections,

    $obs_lookup, $obs_views,
    y
) setup=(x = rand(); y = zeros(n_obs))


inf_view = @view infections[:, 1]
obs_lookup_ind = obs_lookup[1]
n_obs_ind = length(obs_views[1])


@benchmark waning_curve_individual!(
    x, $mu_short, $omega,
    $sigma_long, $sigma_short, $tau,


    $dist_matrix, $time_diff_matrix,
    $(subject_birth_ix[1]),

    $infections,

    $obs_lookup_ind,
    y
) setup=(x = rand(); y = zeros(n_obs_ind))


A = MatrixBernoulli(rand(Uniform(0.0, 1), (10, 10)))

@benchmark Distributions._logpdf($A, M) setup = (M = Matrix{Real}(zeros(Bool, 10, 10)))

M_test = [rand(MatrixBernoulli(rand(Uniform(0.0, 1), (10, 10)))) for i in 1:1000]

@profview_allocs [Distributions._logpdf(A, M_test[i]) for i in 1:1000]


t = Matrix{Real}(zeros(Bool, 10, 10))


data_code = "sim_study_hanam_2018_3"

run_dir = "runs/$(data_code)/"

model_data = load("$run_dir/model_data.hdf5")

obs_df = DataFrame(model_data["observations"])

p = read_model_parameters(model_data)

model = make_waning_model(p, obs_df);

gibbs_sampler = make_gibbs_sampler(model)

sample(model, gibbs_sampler, 2, callback = log_callback);

@profview sample(model, gibbs_sampler, 500);
@profview_allocs sample(model, gibbs_sampler, 1000, callback = log_callback);

@time sample(model, gibbs_sampler, 500, callback = log_callback);

using Cthulhu

@descend model.f(
    model,
    Turing.VarInfo(model),
    Turing.SamplingContext(
        Random.default_rng(), Turing.SampleFromPrior(), Turing.DefaultContext()
    ),
    model.args...,
)




@descend logpdf(a, y)

using BenchmarkTools

a = TitreArrayNormal(rand(300) .+ 3, 0.5, 0.0, 8.0)
y = rand(a)

@benchmark logpdf(a, y) setup = (a = a, y = rand(a))


@profview [logpdf(a, y) for i in 1:10000]



M = MatrixBernoulli(0.15, 100, 100)
M = filldist(Bernoulli(0.15), 100, 100)

y = rand(M)

@benchmark logpdf(M, y) setup = (M = M, y = rand(M))