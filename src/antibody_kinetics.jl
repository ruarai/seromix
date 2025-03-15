
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

    observations,
    obs_lookup, n_obs
)
    y = zeros(typeof(mu_long), n_obs)
    n_t_steps, n_ind = size(infections)
    
    obs_s = view(observations, :, 2)
    
    for i in 1:n_ind
        cumulative_infections = 0.0
        
        for t in 1:n_t_steps
            @inbounds if !infections[t, i]
                continue
            end
            
            cumulative_infections += 1

            # Only process relevant times (after current infection)
            for obs_time in t:n_t_steps
                
                # Skip if no observations match this individual and time
                if !haskey(obs_lookup, (i, obs_time))
                    continue
                end

                matched_indices = obs_lookup[(i, obs_time)]
                
                # Process all matching observations
                for ix_obs in matched_indices
                    @inbounds begin
                        distance = dist_matrix[t, obs_s[ix_obs]]
                        time_diff = obs_time - t
                        
                        # Calculate titre contribution
                        y[ix_obs] += @inline titre_component(
                            mu_long, mu_short, omega,
                            sigma_long, sigma_short, tau,
                            cumulative_infections,
                            distance, time_diff
                        )
                    end
                end
            end
            
        end
    end
    
    return y
end