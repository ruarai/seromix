

function waning_curve!(
    mu_long::T,
    mu_short::T, omega::T,
    sigma_long::T, sigma_short::T,
    tau::T,
    dist_matrix::Matrix{Float64},
    infections::Matrix{Bool},

    obs_lookup::Vector{Dict{Int, Vector{Tuple{Int,Int}}}},
    obs_views::Vector{UnitRange{Int}},
    y::AbstractArray{T}
) where T <: Real
    n_t_steps, n_ind = size(infections)

    
    for i in 1:n_ind
        waning_curve_individual!(
            mu_long, mu_short, omega,
            sigma_long, sigma_short, tau,

            dist_matrix, 

            view(infections, :, i),
            obs_lookup[i],

            n_t_steps,
            view(y, obs_views[i])
        )
    end
    
    return y
end

function waning_curve_individual!(
    mu_long::T,
    mu_short::T, omega::T,
    sigma_long::T, sigma_short::T,
    tau::T,
    dist_matrix::Matrix{Float64},
    infections::AbstractArray{Bool},

    obs_lookup_ind::Dict{Int64, Vector{Tuple{Int64,Int64}}},

    n_t_steps::Int,

    y::AbstractArray{T}
) where T <: Real
    prior_infections = 0.0

    for t in 1:n_t_steps
        @inbounds if !infections[t]
            continue
        end
        
        seniority = max(0.0, 1.0 - tau * prior_infections)
        prior_infections += 1.0
        
        # Process relevant times after infection
        for obs_time in t:n_t_steps
            if !haskey(obs_lookup_ind, (obs_time))
                continue
            end

            matches = obs_lookup_ind[(obs_time)]
            time_diff = obs_time - t
            short_term_time_factor = max(0.0, 1.0 - omega * time_diff)
            
            # Process matching observations
            @inbounds for (s, ix_obs) in matches
                distance = dist_matrix[t, s]
                
                long_term_dist = max(0.0, 1.0 - sigma_long * distance)
                short_term_dist = max(0.0, 1.0 - sigma_short * distance)
                
                # Calculate components
                long_term = mu_long * long_term_dist
                short_term = mu_short * short_term_time_factor * short_term_dist
                
                y[ix_obs] += seniority * (long_term + short_term)
            end
        end
    end
end