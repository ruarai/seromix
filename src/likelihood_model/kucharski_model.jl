
# The Kucharski (2018) model for waning immunity
@model function waning_model_kucharski(
    sp::StaticModelParameters,
    prior_infection_dist::Distribution,
    observed_titre::Vector{Vector{Float64}},
    model_cache::WaningModelCache;

    mixture_importance_sampling::Bool = false,
    use_corrected_titre::Bool = true
)
    mu_long ~ Uniform(0.0, 4.0)
    mu_short ~ Uniform(0.0, 4.0)


    omega ~ Uniform(0.5, 1.0)

    # TODO revert changes
    sigma_long ~ Uniform(0.0, 1.0)
    sigma_short ~ Uniform(0.0, 1.0)

    tau ~ Uniform(0.0, 1.0)

    obs_sd ~ Uniform(1.0, 3.0)

    params = (; mu_long, mu_short, omega, sigma_long, sigma_short, tau, obs_sd)

    infections ~ prior_infection_dist

    @addlogprob! general_waning_likelihood(
        params, infections, observed_titre,

        sp, model_cache, __context__,

        individual_waning_kucharski!;
        use_corrected_titre = use_corrected_titre,
        mixture_importance_sampling = mixture_importance_sampling
    )
end

function individual_waning_kucharski!(
    params,
    infections::AbstractArray{Bool},
    latent_titre::AbstractArray{Float64},
    ix_subject::Int,

    sp::StaticModelParameters,
    model_cache::WaningModelCache
)
    n_t_steps = length(infections)

    prior_infections = 0.0

    for ix_t in max(1, sp.subject_birth_ix[ix_subject]):n_t_steps
        @inbounds if !infections[ix_t]
            continue
        end
        
        seniority = max(0.0, 1.0 - params.tau * prior_infections)
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
                
                latent_titre[ix_obs] += seniority * (long_term + short_term)
            end
        end
    end
end

# Used by pointwise_likelihood()
turing_function_to_waning_function(model_f::typeof(waning_model_kucharski)) = individual_waning_kucharski!
