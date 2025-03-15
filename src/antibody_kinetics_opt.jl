

function waning_curve_optimised(
    mu_long::T,
    mu_short::T, omega::T,
    sigma_long::T, sigma_short::T,
    tau::T,
    dist_matrix::Matrix{Float64},
    infections::Matrix{Bool},
    observations::Matrix{Int},
    obs_lookup::Dict{Tuple{Int,Int},Vector{Int}}, n_obs::Int
) where T <: Real
    y = zeros(T, n_obs)
    n_t_steps, n_ind = size(infections)
    
    # Pre-allocate constants
    one_t = one(T)
    zero_t = zero(T)
    
    # Direct reference to avoid view allocation
    obs_s = observations[:, 2]
    
    for i in 1:n_ind
        cumulative_infections = zero_t
        
        for t in 1:n_t_steps
            @inbounds if !infections[t, i]
                continue
            end
            
            cumulative_infections += one_t
            seniority = max(zero_t, one_t - tau * (cumulative_infections - one_t))
            
            # Process relevant times after infection
            for obs_time in t:n_t_steps
                # Avoid tuple allocation in hot loop
                key = (i, obs_time)
                if !haskey(obs_lookup, key)
                    continue
                end

                matched_indices = obs_lookup[key]
                time_diff = obs_time - t
                short_term_time_factor = max(zero_t, one_t - omega * time_diff)
                
                # Process matching observations
                @inbounds for ix_obs in matched_indices
                    distance = dist_matrix[t, obs_s[ix_obs]]
                    
                    # Pre-compute distance factors
                    long_term_dist = max(zero_t, one_t - sigma_long * distance)
                    short_term_dist = max(zero_t, one_t - sigma_short * distance)
                    
                    # Calculate components
                    long_term = mu_long * long_term_dist
                    short_term = mu_short * short_term_time_factor * short_term_dist
                    
                    y[ix_obs] += seniority * (long_term + short_term)
                end
            end
        end
    end
    
    return y
end