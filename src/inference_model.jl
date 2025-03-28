


struct FixedModelParameters
    n_t_steps::Int
    n_subjects::Int

    antigenic_distances::Matrix{Float64}

    time_diff_matrix::Matrix{Float64}

    subject_birth_ix::Vector{Int}
end



@model function waning_model(
    model_parameters::FixedModelParameters,

    obs_lookup, obs_views,
    n_max_ind_obs::Int,

    observed_titre::Vector{Vector{Float64}}     
)
    # Should be constrained to what we expect y to be 
    # immediately post-infection
    # Otherwise this becomes non-identifiable against getting lots
    # of subsequent infections.
    mu_sum ~ Uniform(3.0, 6.0)

    mu_long ~ Uniform(0.0, 6.0)

    mu_short = mu_sum - mu_long

    # omega ~ Truncated(LogNormal(-1.0, 0.5), 0, 3)
    omega = convert(typeof(mu_long), 0.75)

    sigma_long ~ Uniform(0, 1)
    sigma_short ~ Uniform(0, 1)

    tau ~ Uniform(0, 1)

    # infections = Matrix{Bool}(undef, model_parameters.n_t_steps, model_parameters.n_subjects)

    # for ix_t in 1:model_parameters.n_t_steps, ix_subject in 1:model_parameters.n_subjects
    #     infections[ix_t, ix_subject] ~ Bernoulli(0.1)
    # end

    infections ~ filldist(Bernoulli(0.1), model_parameters.n_t_steps, model_parameters.n_subjects)

    obs_sigma = convert(typeof(mu_long), 1.5)
    obs_min = convert(typeof(mu_long), const_titre_min)
    obs_max = convert(typeof(mu_long), const_titre_max)

    context = DynamicPPL.leafcontext(__context__)
    if context isa IndividualSubsetContext
        subset_context::IndividualSubsetContext = context

        ix_subject = subset_context.ix
        n_obs_subset = length(obs_views[ix_subject])

        y_pred = zeros(typeof(mu_long), n_obs_subset)

        waning_curve_individual!(
            mu_long, mu_short, omega,
            sigma_long, sigma_short, tau,

            model_parameters.antigenic_distances, model_parameters.time_diff_matrix,
            model_parameters.subject_birth_ix[ix_subject],

            AbstractArray{Bool}(view(infections, :, ix_subject)),

            obs_lookup[ix_subject], 
            y_pred
        )

        observed_titre[ix_subject] ~ TitreArrayNormal(y_pred, obs_sigma, obs_min, obs_max)
    else
        y_pred_mem = zeros(typeof(mu_long), n_max_ind_obs)

        for ix_subject in 1:model_parameters.n_subjects
            n_obs_subset = length(obs_views[ix_subject])
            
            y_pred = view(y_pred_mem, 1:n_obs_subset)
            fill!(y_pred, 0.0)

            # y_pred = zeros(typeof(mu_long), n_obs_subset)

            waning_curve_individual!(
                mu_long, mu_short, omega,
                sigma_long, sigma_short, tau,

                model_parameters.antigenic_distances, model_parameters.time_diff_matrix,
                model_parameters.subject_birth_ix[ix_subject],

                AbstractArray{Bool}(view(infections, :, ix_subject)),

                obs_lookup[ix_subject], 
                y_pred
            )

            observed_titre[ix_subject] ~ TitreArrayNormal(y_pred, obs_sigma, obs_min, obs_max)
        end
    end
end