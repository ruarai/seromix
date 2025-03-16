
include("dependencies.jl")

t_steps = 1:4
n_t_steps = length(t_steps)
n_ind = 50

mu_long = 3.0
mu_short = 5.0
omega = 0.05
sigma_long = 0.2
sigma_short = 0.2

tau = 0.3

dist_matrix = [abs(i - j) for i in 1:1.0:n_t_steps, j in 1:1.0:n_t_steps]

infections = zeros(Bool, n_t_steps, n_ind)

infections = rand(Bernoulli(0.5), (n_t_steps, n_ind))

infections_df = DataFrame(stack([[i[1], i[2]] for i in findall(infections)])', :auto)
complete_obs = expand_grid(s = 1:n_t_steps, t = 1:n_t_steps, i = 1:n_ind)

# Ensure dataframe is sorted at this stage
n_obs = nrow(complete_obs)

obs_views = make_obs_views(complete_obs)
obs_lookup = make_obs_lookup_individuals(complete_obs)

y = zeros(n_obs)

waning_curve!(
    mu_long, mu_short, omega,
    sigma_long, sigma_short, tau,

    dist_matrix,
    infections,

    obs_lookup, obs_views,
    y
)

using BenchmarkTools

@benchmark waning_curve!(
    x, $mu_short, $omega,
    $sigma_long, $sigma_short, $tau,

    $dist_matrix,
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