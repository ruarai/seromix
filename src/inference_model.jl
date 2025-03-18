

@model function waning_model(
    dist_matrix, time_diff_matrix,
    subject_birth_ix,
    
    n_ind, n_t_steps,

    obs_lookup, obs_views,
    n_obs,

    observed_titre     
)
    min_inf_titre = 3.0
    max_inf_titre = 6.0

    # Should be constrained to what we expect y to be 
    # immediately post-infection
    # Otherwise this becomes non-identifiable against getting lots
    # of subsequent infections.
    mu_sum ~ Truncated(LogNormal(1.0, 1.0), min_inf_titre, max_inf_titre)

    # Should be such that mu_short > 0 later
    mu_long ~ Truncated(LogNormal(1.0, 1.0), 0, max_inf_titre)

    mu_short = mu_sum - mu_long

    # omega ~ Truncated(LogNormal(-1.0, 0.5), 0, 3)
    omega = convert(typeof(mu_long), 0.8)

    sigma_long ~ Truncated(LogNormal(-2.0, 1.0), 0, 1)
    sigma_short ~ Truncated(LogNormal(-2.0, 1.0), 0, 1)

    tau ~ LogNormal(-2.0, 1.0)

    # Maybe replace with a matrix-variate distribution?
    infections ~ filldist(Bernoulli(0.2), n_t_steps, n_ind)

    context = DynamicPPL.leafcontext(__context__)

    obs_sigma = 1.3

    if context isa IndividualSubsetContext
        subset_context::IndividualSubsetContext = context

        ix_subject = subset_context.ix
        n_obs_subset = length(obs_views[ix_subject])

        y_pred = zeros(typeof(mu_long), n_obs_subset)

        waning_curve_individual!(
            mu_long, mu_short, omega,
            sigma_long, sigma_short, tau,

            dist_matrix, time_diff_matrix,
            subject_birth_ix[ix_subject],
            AbstractArray{Bool}(view(infections, :, ix_subject)),

            obs_lookup[ix_subject], n_t_steps,

            y_pred
        )

        # Replace with one-liner as in below?
        for (i, ix_obs) in enumerate(obs_views[ix_subject])
            observed_titre[ix_obs] ~ Normal(y_pred[i], obs_sigma)
        end
    else
        y_pred = zeros(typeof(mu_long), n_obs)

        waning_curve!(
            mu_long, mu_short, omega,
            sigma_long, sigma_short, tau,
    
            dist_matrix, time_diff_matrix,
            subject_birth_ix,
            Matrix{Bool}(infections), # Is this necessary?
    
            obs_lookup, obs_views,
            y_pred
        )
    
        observed_titre ~ MvNormal(y_pred, obs_sigma * I)
    end
end