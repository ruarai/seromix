
function titre_component(
    mu_long, mu_short, omega,
    sigma_long, sigma_short, tau,
    infection_number,
    distance, t
)
    seniority = max(0, 1 - tau * (infection_number - 1))
    long_term = mu_long * max(0, 1 - sigma_long * distance)
    short_term = mu_short * max(0, 1 - omega * t) * max(0, 1 - sigma_short * distance)
    return seniority * (long_term + short_term)
end

function waning_curve(
    mu_long,
    mu_short, omega,
    sigma_long, sigma_short,
    tau,
    dist_matrix,
    infections,
    observations
)
    n_obs = size(observations, 2)
    y = zeros(typeof(mu_long), n_obs)
    n_t_steps, n_ind = size(infections)
    
    # Pre-extract observation data for faster access
    obs_t = view(observations, 1, :)
    obs_s = view(observations, 2, :)
    obs_i = view(observations, 3, :)
    
    for i in 1:n_ind
        cumulative_infections = 0
        
        for t in 1:n_t_steps
            @inbounds if !infections[t, i]
                continue
            end
            
            cumulative_infections += 1
            
            for ix_obs in 1:n_obs
                @inbounds begin
                    # Skip if observation doesn't match this individual or is before infection
                    if obs_i[ix_obs] != i || obs_t[ix_obs] < t
                        continue
                    end
                    
                    distance = dist_matrix[t, obs_s[ix_obs]]
                    time_diff = obs_t[ix_obs] - t
                    
                    # Calculate titre contribution
                    y[ix_obs] += titre_component(
                        mu_long, mu_short, omega,
                        sigma_long, sigma_short, tau,
                        cumulative_infections,
                        distance, time_diff
                    )
                end
            end
        end
    end
    
    return y
end