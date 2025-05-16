


struct FixedModelParameters
    n_t_steps::Int
    n_subjects::Int

    antigenic_distances::Matrix{Float64}

    time_diff_matrix::Matrix{Float64}

    subject_birth_ix::Vector{Int}
end

function make_waning_model(
    model_parameters::FixedModelParameters,
    obs_df::DataFrame
)
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
    mu_long ~ Uniform(0.0, 10.0)
    mu_short ~ Uniform(0.0, 10.0)


    omega ~ Uniform(0.0, 1.0)

    sigma_long ~ Uniform(0.0, 10.0)
    sigma_short ~ Uniform(0.0, 10.0)

    tau ~ Uniform(0.0, 10.0)

    infections ~ MatrixBernoulli(0.15, model_parameters.n_t_steps, model_parameters.n_subjects)

    obs_sd ~ Uniform(0.0, 10.0)
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