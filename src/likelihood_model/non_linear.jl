@model function waning_model_non_linear(
    model_parameters::FixedModelParameters,
    prior_infection_dist::Distribution,
    observed_titre::Vector{Vector{Float64}},
    model_cache::WaningModelCache;

    mixture_importance_sampling::Bool = false,
    use_corrected_titre::Bool = true
)
    mu_long ~ Uniform(0.0, 10.0)
    mu_long_mult ~ Uniform(0.0, 10.0)

    mu_short ~ Uniform(0.0, 10.0)
    mu_short_mult ~ Uniform(0.0, 10.0)


    omega ~ Uniform(0.0, 1.0)
    omega_mult ~ Uniform(0.0, 1.0)

    sigma_long ~ Uniform(0.0, 10.0)
    sigma_long_mult ~ Uniform(0.0, 10.0)

    sigma_short ~ Uniform(0.0, 10.0)
    sigma_short_mult ~ Uniform(0.0, 10.0)

    tau ~ Uniform(0.0, 1.0)
    beta ~ Uniform(0.0, 1.0)

    obs_sd ~ Uniform(1.0, 10.0)

    intercept ~ Uniform(-8.0, 1.0)

    params = (; 
        mu_long, mu_long_mult, mu_short, mu_short_mult,
        omega, omega_mult, sigma_long, sigma_long_mult,
        sigma_short, sigma_short_mult, beta, tau, obs_sd, intercept
    )

    infections ~ prior_infection_dist

    log_likelihood = general_waning_likelihood(
        params,
        infections,
        model_parameters,
        observed_titre,
        model_cache,
        individual_waning_non_linear!,
        DynamicPPL.leafcontext(__context__);
        use_corrected_titre = use_corrected_titre,
        mixture_importance_sampling = mixture_importance_sampling
    )

    @addlogprob! log_likelihood
end

function individual_waning_non_linear!(
    params,
    dist_matrix::Matrix{Float64},
    time_diff_matrix::Matrix{Float64},
    subject_birth_ix::Int,
    infections::AbstractArray{Bool},

    obs_lookup_strain,
    obs_lookup_ix,

    y::AbstractArray{Float64}
)
    n_t_steps = length(infections)

    prior_infections = 0.0

    y .= 1.0

    for ix_t in max(1, subject_birth_ix):n_t_steps
        @inbounds if !infections[ix_t]
            continue
        end
        
        age = ix_t - subject_birth_ix
        
        age_effect = max(0.0, 1.0 - age * params.beta)
        seniority = max(0.0, 1.0 - params.tau * prior_infections)

        prior_infections += 1.0

        # Process relevant times after infection
        for ix_t_obs in ix_t:n_t_steps
            if !haskey(obs_lookup_strain, ix_t_obs)
                continue
            end

            matches_strain = obs_lookup_strain[ix_t_obs]
            matches_ix = obs_lookup_ix[ix_t_obs]

            time_diff = time_diff_matrix[ix_t_obs, ix_t]
            # short_term_time_factor = 1.0 - params.omega * time_diff 
            
            @turbo for i in eachindex(matches_strain)
                ix_obs_strain = matches_strain[i]
                ix_obs = matches_ix[i]

                distance = dist_matrix[ix_t, ix_obs_strain]
                
                long_term_dist = 1.0 - params.sigma_long * distance
                # short_term_dist = 1.0 - params.sigma_short * distance
                
                long_term = params.mu_long * long_term_dist
                short_term = params.mu_short * (1 - params.omega * time_diff - params.sigma_short * distance)

                y[ix_obs] += 2 ^ (seniority * age_effect * long_term) +
                    2 ^ (seniority * age_effect * short_term)
            end
        end
    end

    # Return to log-scale and add intercept
    y .= log2.(y) .+ params.intercept

    # Repeat loop for short-term effect (multiplicative, so additive on log-scale)
    prior_infections = 0.0

    for ix_t in max(1, subject_birth_ix):n_t_steps
        @inbounds if !infections[ix_t]
            continue
        end
        
        age = ix_t - subject_birth_ix
        
        age_effect = max(0.0, 1.0 - age * params.beta)
        seniority = max(0.0, 1.0 - params.tau * prior_infections)

        prior_infections += 1.0

        # Process relevant times after infection
        for ix_t_obs in ix_t:n_t_steps
            if !haskey(obs_lookup_strain, ix_t_obs)
                continue
            end

            matches_strain = obs_lookup_strain[ix_t_obs]
            matches_ix = obs_lookup_ix[ix_t_obs]

            time_diff = time_diff_matrix[ix_t_obs, ix_t]
            short_term_time_factor = max(0.0, 1.0 - params.omega_mult * time_diff)
            
            @turbo for i in eachindex(matches_strain)
                ix_obs_strain = matches_strain[i]
                ix_obs = matches_ix[i]

                distance = dist_matrix[ix_t, ix_obs_strain]
                
                long_term_dist = max(0.0, 1.0 - params.sigma_long_mult * distance)
                short_term_dist = max(0.0, 1.0 - params.sigma_short_mult * distance)

                long_term = params.mu_long_mult * long_term_dist
                short_term = params.mu_short_mult * short_term_time_factor * short_term_dist

                y[ix_obs] += seniority * age_effect * (long_term + short_term)
            end
        end
    end
end

# Used by pointwise_likelihood()
turing_function_to_waning_function(model_f::typeof(waning_model_non_linear)) = individual_waning_non_linear!