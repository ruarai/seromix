
function waning_curve_individual_linear!(
    mu_add::T, mu_mult::T, dist_scale_add::T, dist_scale_mult::T,
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
        
        # Process relevant times after infection
        for ix_t_obs in ix_t:n_t_steps
            if !haskey(obs_lookup_ind, ix_t_obs)
                continue
            end

            matches = obs_lookup_ind[ix_t_obs]
            time_diff = time_diff_matrix[ix_t_obs, ix_t]
            
            @inbounds for (ix_obs_strain, ix_obs) in matches
                dist_add = dist_matrix[ix_t, ix_obs_strain] * dist_scale_add
                dist_mult = dist_matrix[ix_t, ix_obs_strain] * dist_scale_mult
                
                additive_effect = mu_add - dist_add
                multiplicative_effect = mu_mult * max(0.0, 1  - dist_mult)

                if prior_infections > 0
                    y[ix_obs] *= 2 ^ multiplicative_effect
                end
                
                y[ix_obs] += 2 ^ additive_effect
            end
        end

        prior_infections += 1.0
    end
end


