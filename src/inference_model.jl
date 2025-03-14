

@model function waning_model(
    dist_matrix,
    
    n_ind, n_t_steps,
    
    obs_df, obs_titre
)
    mu_long ~ LogNormal(1.0, 1.0)

    # Should be constrained to what we expect y to be 
    # immediately post-infection
    # Otherwise this becomes non-identifiable against getting lots
    # of subsequent infections.
    mu_sum ~ Truncated(LogNormal(1.0, 1.0), 6.5, 10)

    mu_short = mu_sum - mu_long

    # omega ~ LogNormal(-3.0, 1.0)
    omega = 0.05

    sigma_long ~ Truncated(LogNormal(-2.0, 1.0), 0, 1)
    sigma_short ~ Truncated(LogNormal(-2.0, 1.0), 0, 1)

    # Also has non-identifiability problems with the
    # infection matrix !
    tau ~ LogNormal(-2.0, 1.0)

    infections ~ filldist(Bernoulli(0.2), n_t_steps, n_ind)

    y_pred = waning_curve(
        mu_long, mu_short, omega,
        sigma_long, sigma_short, tau,

        dist_matrix,
        infections,
        obs_df
    )

    obs_titre ~ MvNormal(y_pred, 0.3 * I)
end