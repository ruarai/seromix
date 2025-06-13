
# The Kucharski (2018) model for waning immunity
@model function waning_model_kucharski(
    model_parameters::FixedModelParameters,
    prior_infection_dist::Distribution,
    use_corrected_titre::Bool,

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

    obs_sd ~ Uniform(0.0, 10.0)
    obs_min = convert(typeof(mu_long), const_titre_min)
    obs_max = convert(typeof(mu_long), const_titre_max)

    infections ~ prior_infection_dist

    # If we're in an "IndividualSubsetContext", only calculate the
    # likelihood over a single subject
    # Otherwise, calculate across all subjects
    context = DynamicPPL.leafcontext(__context__)
    subjects_to_process = if context isa IndividualSubsetContext
        SA[context.ix]
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
            AbstractArray{Bool}(view(infections, :, ix_subject)),
            obs_lookup[ix_subject], 
            y_pred
        )

         observed_titre[ix_subject] ~ TitreArrayNormal(y_pred, obs_sd, obs_min, obs_max, use_corrected_titre)
    end
end


function ppd_kucharski(chain, model_parameters, n_ppd_subjects, n_draws)
    # TODO can this be refactored over an arbitrary vector of ix_subject?
    obs_df = expand_grid(
        ix_t_obs = 1:model_parameters.n_t_steps, 
        ix_strain = 1:model_parameters.n_t_steps, 
        ix_subject = 1:n_ppd_subjects,
        observed_titre = 0.0,
        ix_draw = 1:n_draws,
    )

    obs_df_grouped = groupby(obs_df, :ix_draw)


    col_names = names(chain)
    ix_infections = findall(s -> startswith(s, "infections"), col_names)

    for ix_draw in 1:n_draws
        ix_sample = sample(1:nrow(chain)) 

        draw = chain[ix_sample,:]
        
        # TODO cache across ix_draw
        obs_lookup = make_obs_lookup(obs_df_grouped[ix_draw])
        obs_view = make_obs_views(obs_df_grouped[ix_draw])

        infections = reshape(Vector(chain[ix_sample, ix_infections]), model_parameters.n_t_steps, model_parameters.n_subjects)

        for ix_subject in 1:n_ppd_subjects
            waning_curve_individual!(
                draw.mu_long, draw.mu_short, draw.omega,
                draw.sigma_long, draw.sigma_short, draw.tau,

                model_parameters.antigenic_distances, model_parameters.time_diff_matrix,
                model_parameters.subject_birth_ix[ix_subject],

                AbstractArray{Bool}(view(infections, :, ix_subject)),

                obs_lookup[ix_subject], 
                view(obs_df_grouped[ix_draw].observed_titre, obs_view[ix_subject])
            )
        end
    end

    return DataFrame(obs_df_grouped)
end