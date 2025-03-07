

@model function waning_model(
    t_obs, 
    
    t_infected,
    
    titer_obs
)
    mu ~ Normal(3.0, 1.0)
    omega ~ LogNormal(-3.0, 1.0)

    y_pred = waning_curve(mu, omega, t_infected, t_obs)

    titer_obs ~ MvNormal(y_pred, 0.5 * I)
end