
function waning_curve_individual_linear!(
    mu_long::T, mu_short::T, dist_scale::T, r_waning::T,
    dist_matrix::Matrix{Float64},
    time_diff_matrix::Matrix{Float64},
    subject_birth_ix::Int,
    infections::AbstractArray{Bool},

    obs_lookup_ind::Dict{Int64, Vector{Tuple{Int64,Int64}}},

    y::AbstractArray{T}
) where T <: Real
    n_t_steps = length(infections)

    for ix_t in max(1, subject_birth_ix):n_t_steps
        @inbounds if !infections[ix_t]
            continue
        end
        
        # Process relevant times after infection
        for ix_t_obs in ix_t:n_t_steps
            if !haskey(obs_lookup_ind, ix_t_obs)
                continue
            end

            matches = obs_lookup_ind[ix_t_obs]
            time_diff = time_diff_matrix[ix_t_obs, ix_t]
            
            @inbounds for (ix_obs_strain, ix_obs) in matches
                distance_scaled = dist_matrix[ix_t, ix_obs_strain] * dist_scale

                long_term_effect = 2 ^ (mu_long - distance_scaled)
                short_term_effect = 2 ^ (mu_short - distance_scaled - r_waning * time_diff)
                
                y[ix_obs] += 2 ^ long_term_effect + 2 ^ (short_term_effect)
            end
        end
    end
end