

function waning_curve!(
    mu_long::T,
    mu_short::T, omega::T,
    sigma_long::T, sigma_short::T,
    tau::T,
    dist_matrix::Matrix{Float64},
    time_diff_matrix::Matrix{Float64},
    subject_birth_ix::Vector{Int},
    infections::Matrix{Bool},

    obs_lookup::Vector{Dict{Int, Vector{Tuple{Int,Int}}}},
    obs_views::Vector{UnitRange{Int}},
    y::AbstractArray{T}
) where T <: Real
    n_t_steps, n_ind = size(infections)

    
    for ix_subject in 1:n_ind
        waning_curve_individual!(
            mu_long, mu_short, omega,
            sigma_long, sigma_short, tau,

            dist_matrix, time_diff_matrix,
            subject_birth_ix[ix_subject],

            view(infections, :, ix_subject),
            obs_lookup[ix_subject],

            view(y, obs_views[ix_subject])
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
    time_diff_matrix::Matrix{Float64},
    subject_birth_ix::Int,
    infections::AbstractArray{Bool},

    obs_lookup_ind::Dict{Int64, Vector{Tuple{Int64,Int64}}},

    y::AbstractArray{T}
) where T <: Real
    n_t_steps = length(infections)

    prior_infections = 0.0

    for ix_t in max(1, subject_birth_ix):n_t_steps
        @inbounds if !infections[ix_t]
            continue
        end
        
        seniority = max(0.0, 1.0 - tau * prior_infections)
        prior_infections += 1.0
        
        # Process relevant times after infection
        for ix_t_obs in ix_t:n_t_steps
            if !haskey(obs_lookup_ind, ix_t_obs)
                continue
            end

            matches = obs_lookup_ind[ix_t_obs]
            time_diff = time_diff_matrix[ix_t_obs, ix_t]
            short_term_time_factor = max(0.0, 1.0 - omega * time_diff)
            
            @inbounds for (ix_obs_strain, ix_obs) in matches
                distance = dist_matrix[ix_t, ix_obs_strain]
                
                long_term_dist = max(0.0, 1.0 - sigma_long * distance)
                short_term_dist = max(0.0, 1.0 - sigma_short * distance)
                
                long_term = mu_long * long_term_dist
                short_term = mu_short * short_term_time_factor * short_term_dist
                
                y[ix_obs] += seniority * (long_term + short_term)
            end
        end
    end
end