@model function waning_model_age_effect(
    sp::StaticModelParameters,
    prior_infection_dist::Distribution,
    observed_titre::Vector{Vector{Float64}},
    model_cache::WaningModelCache;

    mixture_importance_sampling::Bool = false,
    use_corrected_titre::Bool = true
)
    mu_long ~ Uniform(0.0, 10.0)
    mu_short ~ Uniform(0.0, 10.0)


    omega ~ Uniform(0.0, 1.0)

    sigma_long ~ Uniform(0.0, 1.0)
    sigma_short ~ Uniform(0.0, 1.0)

    beta ~ Uniform(0.0, 1.0)
    tau ~ Uniform(0.0, 1.0)

    tau_cutoff ~ Uniform(0.0, 1.0)

    intercept ~ Uniform(-5.0, 1.0)

    obs_sd ~ Uniform(1.0, 10.0)

    params = (; mu_long, mu_short, omega, sigma_long, sigma_short, beta, tau, tau_cutoff, intercept, obs_sd)

    infections ~ prior_infection_dist

    @addlogprob! general_waning_likelihood(
        params, infections, observed_titre,

        sp, model_cache, __context__,

        individual_waning_age_effect!;
        use_corrected_titre = use_corrected_titre,
        mixture_importance_sampling = mixture_importance_sampling
    )
end


function individual_waning_age_effect!(
    params,
    infections::AbstractArray{Bool},
    latent_titre::AbstractArray{Float64},
    ix_subject::Int,

    sp::StaticModelParameters,
    model_cache::WaningModelCache
)
    n_t_steps = length(infections)
    
    latent_titre .+= params.intercept
    prior_infections = 0.0

    subject_birth_ix = sp.subject_birth_ix[ix_subject]

    for ix_t in max(1, subject_birth_ix):n_t_steps
        @inbounds if !infections[ix_t]
            continue
        end
        
        age = ix_t - subject_birth_ix
        
        age_effect = max(0.0, 1.0 - age * params.beta)

        # Increase 10 here?
        seniority = max(1 - params.tau * (params.tau_cutoff * 15), 1.0 - params.tau * prior_infections)

        prior_infections += 1.0

        # Process relevant times after infection
        for ix_t_obs in ix_t:n_t_steps
            if !haskey(model_cache.obs_lookup_strain[ix_subject], ix_t_obs)
                continue
            end

            matches_strain = model_cache.obs_lookup_strain[ix_subject][ix_t_obs]
            matches_ix = model_cache.obs_lookup_ix[ix_subject][ix_t_obs]

            time_diff = sp.time_diff_matrix[ix_t_obs, ix_t]
            short_term_time_factor = max(0.0, 1.0 - params.omega * time_diff)
            
            @turbo for i in eachindex(matches_strain)
                ix_obs_strain = matches_strain[i]
                ix_obs = matches_ix[i]

                distance = sp.antigenic_distances[ix_t, ix_obs_strain]
                
                long_term_dist = max(0.0, 1.0 - params.sigma_long * distance)
                short_term_dist = max(0.0, 1.0 - params.sigma_short * distance)
                
                long_term = params.mu_long * long_term_dist
                short_term = params.mu_short * short_term_time_factor * short_term_dist
                
                latent_titre[ix_obs] += seniority * age_effect * (long_term + short_term)
            end
        end
    end
end

# Used by pointwise_likelihood()
turing_function_to_waning_function(model_f::typeof(waning_model_age_effect)) = individual_waning_age_effect!