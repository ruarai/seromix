
function make_waning_model(
    model_parameters::FixedModelParameters,
    obs_df::DataFrame
)
    all(diff(obs_df.ix_subject) .>= 0) || throw(ArgumentError("ix_subject in obs_df must be sorted in ascending order."))

    n_max_ind_obs = maximum(length.(make_obs_views(obs_df)))

    individual_titre_obs = [obs_df.observed_titre[v] for v in make_obs_views(obs_df)]
    return waning_model(
        model_parameters,

        make_obs_lookup(obs_df),
        make_obs_views(obs_df),
        n_max_ind_obs,
        individual_titre_obs
    );
end

@model function waning_model(
    model_parameters::FixedModelParameters,

    obs_lookup, obs_views,
    n_max_ind_obs::Int,

    observed_titre::Vector{Vector{Float64}}     
)
    latent_means = [log(2), log(2), 0.0, log(0.5), log(0.5), log(0.1), log(1)]
    latent_basic_params ~ MvNormal(latent_means, I)

    omega_link_dist = Uniform(0, 1)

    mu_long := exp(latent_basic_params[1])
    mu_short := exp(latent_basic_params[2])
    omega := invlink(omega_link_dist, latent_basic_params[3])
    sigma_long := exp(latent_basic_params[4])
    sigma_short := exp(latent_basic_params[5])
    tau := exp(latent_basic_params[6])
    obs_sd := exp(latent_basic_params[7])
    
    obs_min = convert(typeof(mu_long), const_titre_min)
    obs_max = convert(typeof(mu_long), const_titre_max)

    time_effect ~ MvNormal(fill(-1.0, model_parameters.n_t_steps), I)
    subject_effect ~ MvNormal(fill(0.0, model_parameters.n_subjects), I)

    infections ~ MatrixHierarchicalBernoulli(
        time_effect,
        subject_effect,
        model_parameters.n_t_steps,
        model_parameters.n_subjects
    )
    
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

            # TODO remove abstractarray here?
            AbstractArray{Bool}(view(infections, :, ix_subject)),

            obs_lookup[ix_subject], 
            y_pred
        )

        observed_titre[ix_subject] ~ TitreArrayNormal(y_pred, obs_sd, obs_min, obs_max)
    else
        y_pred_mem = zeros(typeof(mu_long), n_max_ind_obs)

        for ix_subject in 1:model_parameters.n_subjects
            n_obs_subset = length(obs_views[ix_subject])
            
            y_pred = view(y_pred_mem, 1:n_obs_subset)
            fill!(y_pred, 0.0)

            waning_curve_individual!(
                mu_long, mu_short, omega,
                sigma_long, sigma_short, tau,

                model_parameters.antigenic_distances, model_parameters.time_diff_matrix,
                model_parameters.subject_birth_ix[ix_subject],

                AbstractArray{Bool}(view(infections, :, ix_subject)),

                obs_lookup[ix_subject], 
                y_pred
            )

            observed_titre[ix_subject] ~ TitreArrayNormal(y_pred, obs_sd, obs_min, obs_max)
        end
    end
end