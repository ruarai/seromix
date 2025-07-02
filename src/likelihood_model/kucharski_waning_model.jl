
# The Kucharski (2018) model for waning immunity
@model function waning_model_kucharski(
    model_parameters::FixedModelParameters,
    prior_infection_dist::Distribution,
    use_corrected_titre::Bool,

    obs_lookup_strain, obs_lookup_ix, obs_views,
    n_max_ind_obs::Int,

    observed_titre::Vector{Vector{Float64}}     
)
    mu_long ~ Uniform(0.0, 10.0)
    mu_short ~ Uniform(0.0, 10.0)


    omega ~ Uniform(0.0, 1.0)

    sigma_long ~ Uniform(0.0, 10.0)
    sigma_short ~ Uniform(0.0, 10.0)

    tau ~ Uniform(0.0, 10.0)

    obs_sd ~ Uniform(0.0, 10.0)
    obs_min = convert(typeof(mu_long), const_titre_min)
    obs_max = convert(typeof(mu_long), const_titre_max)

    infections ~ prior_infection_dist

    # If we're in an "IndividualSubsetContext", only calculate the
    # likelihood over a single subject
    # Otherwise, calculate across all subjects
    context = DynamicPPL.leafcontext(__context__)
    subjects_to_process = if context isa IndividualSubsetContext
        SA[context.ix] # StaticArray with one element (context.ix)
    else
        1:model_parameters.n_subjects
    end

    # Reduce the total memory allocation across the observed titre
    # by calculating across a single pre-allocated array.
    y_pred_buffer = zeros(typeof(mu_long), n_max_ind_obs)

    for ix_subject in subjects_to_process
        n_obs_subject = length(obs_views[ix_subject])
        y_pred = view(y_pred_buffer, 1:n_obs_subject)
        fill!(y_pred, 0.0)

        waning_curve_individual!(
            mu_long, mu_short, omega,
            sigma_long, sigma_short, tau,
            model_parameters.antigenic_distances, model_parameters.time_diff_matrix,
            model_parameters.subject_birth_ix[ix_subject],
            view(infections, :, ix_subject),
            obs_lookup_strain[ix_subject],
            obs_lookup_ix[ix_subject],
            y_pred
        )

         observed_titre[ix_subject] ~ TitreArrayNormal(y_pred, obs_sd, obs_min, obs_max, use_corrected_titre)
    end
end

