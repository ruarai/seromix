

function waning_curve!(
    mu_long::T,
    mu_short::T, omega::T,
    sigma_long::T, sigma_short::T,
    tau::T,
    dist_matrix::Matrix{Float64},
    time_diff_matrix::Matrix{Float64},
    subject_birth_ix::Vector{Int},
    infections::Matrix{Bool},
    obs_lookup_strain,
    obs_lookup_ix,
    obs_views::Vector{UnitRange{Int}},
    y::AbstractArray{T}
) where T <: Real
    for ix_subject in axes(infections, 2)
        waning_curve_individual!(
            mu_long, mu_short, omega,
            sigma_long, sigma_short, tau,

            dist_matrix, time_diff_matrix,
            subject_birth_ix[ix_subject],

            view(infections, :, ix_subject),
            obs_lookup_strain[ix_subject],
            obs_lookup_ix[ix_subject],

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

    obs_lookup_strain,
    obs_lookup_ix,

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
            if !haskey(obs_lookup_strain, ix_t_obs)
                continue
            end

            matches_strain = obs_lookup_strain[ix_t_obs]
            matches_ix = obs_lookup_ix[ix_t_obs]

            time_diff = time_diff_matrix[ix_t_obs, ix_t]
            short_term_time_factor = max(0.0, 1.0 - omega * time_diff)
            
            @turbo for i in eachindex(matches_strain)
                ix_obs_strain = matches_strain[i]
                ix_obs = matches_ix[i]

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


# Necessary to get logjoint to work
function waning_curve_individual!(
    mu_long::T,
    mu_short::T, omega::T,
    sigma_long::T, sigma_short::T,
    tau::T,
    dist_matrix::Matrix{Float64},
    time_diff_matrix::Matrix{Float64},
    subject_birth_ix::Int,
    infections::AbstractArray{Float64},

    obs_lookup_strain,
    obs_lookup_ix,

    y::AbstractArray{T}
) where T <: Real
    waning_curve_individual!(
        mu_long, mu_short, omega, sigma_long, sigma_short,
        tau, dist_matrix, time_diff_matrix, subject_birth_ix,
        convert.(Bool, infections),
        obs_lookup_strain, obs_lookup_ix,
        y
    )
end