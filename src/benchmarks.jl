
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

infections = zeros(Bool, n_t_steps, n_subjects)

infections = rand(Bernoulli(0.5), (n_t_steps, n_subjects))

infections_df = DataFrame(stack([[i[1], i[2]] for i in findall(infections)])', :auto)
complete_obs = expand_grid(ix_strain = 1:n_t_steps, ix_t_obs = 1:n_t_steps, ix_subject = 1:n_subjects)

# Ensure dataframe is sorted at this stage
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

@benchmark waning_curve!(
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

    $dist_matrix, $inf_view,

    $obs_lookup_ind, $n_t_steps,
    y
) setup=(x = rand(); y = zeros(n_obs_ind))


A = MatrixBernoulli(rand(Uniform(0.0, 1), (10, 10)))

@benchmark Distributions._logpdf($A, M) setup = (M = Matrix{Real}(zeros(Bool, 10, 10)))

M_test = [rand(MatrixBernoulli(rand(Uniform(0.0, 1), (10, 10)))) for i in 1:1000]

@profview_allocs [Distributions._logpdf(A, M_test[i]) for i in 1:1000]


t = Matrix{Real}(zeros(Bool, 10, 10))