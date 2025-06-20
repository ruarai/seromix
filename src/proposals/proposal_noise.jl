
function proposal_noise(
    rng, theta::AbstractVector{Bool}, ix_subject::Int, n_t_steps::Int
)
    ix_start = (ix_subject - 1) * n_t_steps + 1
    ix_end = ix_start + n_t_steps - 1
    theta_view = @view theta[ix_start:ix_end]

    n_swap = clamp(rand(rng, Poisson(10.0)), 1, length(theta_view))
    ix_swap = sample(rng, eachindex(theta_view), n_swap, replace = false)

    log_hastings_ratio = 0.0

    return ix_swap, log_hastings_ratio
end