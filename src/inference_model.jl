

@model function waning_model(
    dist_matrix,
    
    n_ind, n_t_steps,

    obs_lookup, obs_views,
    n_obs,

    obs_titre     
)
    min_titer_inf = 3.0
    max_titer_inf = 6.0

    # Should be constrained to what we expect y to be 
    # immediately post-infection
    # Otherwise this becomes non-identifiable against getting lots
    # of subsequent infections.
    mu_sum ~ Truncated(LogNormal(1.0, 1.0), min_titer_inf, max_titer_inf)

    # Should be such that mu_short > 0 later
    mu_long ~ Truncated(LogNormal(1.0, 1.0), 0, max_titer_inf)

    mu_short = mu_sum - mu_long

    # omega ~ LogNormal(-3.0, 1.0)
    omega = convert(typeof(mu_long), 0.8)

    sigma_long ~ Truncated(LogNormal(-2.0, 1.0), 0, 1)
    sigma_short ~ Truncated(LogNormal(-2.0, 1.0), 0, 1)

    # Also has non-identifiability problems with the
    # infection matrix !
    tau ~ LogNormal(-2.0, 1.0)

    infections ~ filldist(Bernoulli(0.2), n_t_steps, n_ind)

    context = DynamicPPL.leafcontext(__context__)

    if context isa IndividualSubsetContext
        subset_context::IndividualSubsetContext = context

        ix_ind = subset_context.ix
        n_obs_subset = length(obs_views[ix_ind])

        y_pred = zeros(typeof(mu_long), n_obs_subset)

        waning_curve_individual!(
            mu_long, mu_short, omega,
            sigma_long, sigma_short, tau,

            dist_matrix, 
            AbstractArray{Bool}(view(infections, :, ix_ind)),

            obs_lookup[ix_ind], n_t_steps,

            y_pred
        )

        for (i, ix_obs) in enumerate(obs_views[ix_ind])
            obs_titre[ix_obs] ~ Normal(y_pred[i], 0.5)
        end

    else
        y_pred = zeros(typeof(mu_long), n_obs)

        waning_curve!(
            mu_long, mu_short, omega,
            sigma_long, sigma_short, tau,
    
            dist_matrix,
            Matrix{Bool}(infections), # Is this necessary?
    
            obs_lookup, obs_views,
            y_pred
        )
    
        obs_titre ~ MvNormal(y_pred, 0.5 * I)
    end

end