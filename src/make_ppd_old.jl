
## Terrible function!
function make_ppd(
    chain, n_samples,

    params
)
    ppd_subjects = 1:params.n_subjects

    draw_dfs = Array{DataFrame}(undef, n_samples)

    for ix_draw in 1:n_samples
        chain_sample = sample(chain, 1)

        draw_dfs[ix_draw] = expand_grid(
            ix_t_obs = 1:params.n_t_steps, 
            ix_strain = 1:params.n_t_steps, 
            ix_subject = ppd_subjects,
            draw = ix_draw,
            observed_titre = 0.0
        )

        obs_lookup = make_obs_lookup(draw_dfs[ix_draw])
        obs_views = make_obs_views(draw_dfs[ix_draw])


        n_obs_subset = length(obs_views[1])

        p = get_params(chain_sample)
        
        pred_inf = ntuple_to_matrix(p.infections, params.n_t_steps, params.n_subjects)

        for ix_subject in ppd_subjects

            y_pred = zeros(n_obs_subset)

            waning_curve_individual!(
                p.mu_long[1], p.mu_sum[1] - p.mu_long[1], 0.8,
                p.sigma_long[1], p.sigma_short[1], p.tau[1],

                params.antigenic_distances, params.time_diff_matrix,
                params.subject_birth_ix[ix_subject],
                AbstractArray{Bool}(view(pred_inf, :, ix_subject)),

                obs_lookup[ix_subject],

                y_pred
            )

            for (i, ix_obs) in enumerate(obs_views[ix_subject])
                draw_dfs[ix_draw].observed_titre[ix_obs] = y_pred[i]
            end
        end
    end

    ppd_obs = vcat(draw_dfs...)
    
    return ppd_obs
end