

@model function waning_model(
    dist_matrix,
    infections,
    n_strain,
    
    obs_df, obs_titre
)
    mu_long ~ Normal(3.0, 1.0)
    mu_short ~ Normal(3.0, 1.0)

    omega ~ LogNormal(-3.0, 1.0)

    sigma_short ~ LogNormal(-2.0, 1.0)
    sigma_long ~ LogNormal(-2.0, 1.0)

    tau ~ LogNormal(-2.0, 1.0)

    unique_t_steps = unique(obs_df.t)

    y_pred = waning_curve(
        mu_long,
        mu_short, omega,
        sigma_long, sigma_short,
        tau,
    
    
        dist_matrix,
        infections, n_strain, 
        unique_t_steps
    )

    y_pred_array = y_pred[CartesianIndex.(zip(obs_df.ix_strain, obs_df.ix_t))]

    obs_titre ~ MvNormal(y_pred_array, 0.5)
end