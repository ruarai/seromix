


include("dependencies.jl")

t_steps = 1:4
n_t_steps = length(t_steps)
n_ind = 500

mu_long = 3.0
mu_short = 5.0
omega = 0.05
sigma_long = 0.2
sigma_short = 0.2

tau = 0.3

dist_matrix = [abs(i - j) for i in 1:1.0:n_t_steps, j in 1:1.0:n_t_steps]

infections = zeros(Bool, n_t_steps, n_ind)

infections = rand(Bernoulli(0.2), (n_t_steps, n_ind))

infections_df = DataFrame(stack([[i[1], i[2]] for i in findall(infections)])', :auto)
write_parquet("data/infections_df.parquet", infections_df)

complete_obs = expand_grid(t = 1:n_t_steps, s = 1:n_t_steps, i = 1:n_ind)

complete_obs.y .= waning_curve_optimised(
    mu_long, mu_short, omega,
    sigma_long, sigma_short, tau,

    dist_matrix,
    infections,

    make_obs_matrix(complete_obs),
    make_obs_lookup(complete_obs), nrow(complete_obs)
)


using BenchmarkTools


n_obs = nrow(complete_obs)

obs_lookup = make_obs_lookup(complete_obs)
obs_matrix = make_obs_matrix(complete_obs)

@time [waning_curve_optimised(
    i, mu_short, omega,
    sigma_long, sigma_short, tau,

    dist_matrix,
    infections,
    obs_matrix,
    obs_lookup, n_obs
) for i in 1:1.0:100];


@benchmark waning_curve_optimised(
    x, $mu_short, $omega,
    $sigma_long, $sigma_short, $tau,

    $dist_matrix,
    $infections,

    $obs_matrix,
    $obs_lookup, $n_obs
) setup=(x = rand())