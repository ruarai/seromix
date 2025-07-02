
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
        # TODO fix with new lookup scheme
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



function pointwise_likelihood_kucharski(chain, model_data)
    model_data = convert_model_data(model_data)
    p = read_model_parameters(model_data)

    obs_df = DataFrame(model_data["observations"])

    col_names = names(chain)
    ix_infections = findall(s -> startswith(s, "infections"), col_names)

    obs_lookup_strain, obs_lookup_ix = make_obs_lookup(obs_df)
    obs_views = make_obs_views(obs_df)

    logp = zeros(nrow(obs_df), nrow(chain))

    y_pred_buffer = zeros(maximum([length(v) for v in obs_views]))

    @showprogress for ix_sample in 1:nrow(chain)
        draw = chain[ix_sample, :]
        

        infections = reshape(Vector(chain[ix_sample, ix_infections]), p.n_t_steps, p.n_subjects)

        for ix_subject in 1:p.n_subjects
            n_obs_subject = length(obs_views[ix_subject])
            y_pred = view(y_pred_buffer, 1:n_obs_subject)
            fill!(y_pred, 0.0)

            waning_curve_individual!(
                draw.mu_long, draw.mu_short, draw.omega,
                draw.sigma_long, draw.sigma_short, draw.tau,

                p.antigenic_distances, p.time_diff_matrix,
                p.subject_birth_ix[ix_subject],

                AbstractArray{Bool}(view(infections, :, ix_subject)),

                obs_lookup_strain[ix_subject], obs_lookup_ix[ix_subject], 
                y_pred
            )

            @inbounds for ix_obs in 1:n_obs_subject
                ix_obs_absolute = obs_views[ix_subject][ix_obs]
                logp[ix_obs_absolute, ix_sample] = titre_logpdf_component(
                    obs_df.observed_titre[ix_obs_absolute],
                    y_pred[ix_obs],
                    draw.obs_sd,
                    const_titre_min,
                    const_titre_max
                )
            end
        end
    end

    return logp
end