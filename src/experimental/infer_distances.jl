

@model function waning_model_free(
    model_parameters::FixedModelParameters,
    prior_infection_dist::Distribution,

    obs_lookup, obs_views,
    n_max_ind_obs::Int,

    observed_titre::Vector{Vector{Float64}}     
)
    latent_means = [log(2), log(2), 0.0, log(0.1), log(1)]
    latent_basic_params ~ MvNormal(latent_means, I)

    omega_link_dist = Uniform(0, 1)

    mu_long := exp(latent_basic_params[1])
    mu_short := exp(latent_basic_params[2])
    omega := invlink(omega_link_dist, latent_basic_params[3])
    sigma_long := convert(typeof(mu_long), 0.15)
    sigma_short := convert(typeof(mu_long), 0.05)
    tau := exp(latent_basic_params[4])
    obs_sd := exp(latent_basic_params[5])


    obs_min = convert(typeof(mu_long), const_titre_min)
    obs_max = convert(typeof(mu_long), const_titre_max)

    strain_gaps ~ MvNormal(
        zeros(typeof(mu_long), model_parameters.n_t_steps - 1),
        create_special_matrix_iterative(model_parameters.n_t_steps - 1)
    )

    strain_locations = vcat([0], cumsum(exp.(strain_gaps)))
    
    antigenic_distances = zeros(typeof(mu_long), model_parameters.n_t_steps, model_parameters.n_t_steps)

    for j in 1:model_parameters.n_t_steps
        for i in 1:j
            loc_i = strain_locations[i]
            loc_j = strain_locations[j]

            antigenic_distances[i, j] = sqrt((loc_i - loc_j)^2)
            antigenic_distances[j, i] = antigenic_distances[i, j]
        end
    end

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
            antigenic_distances, model_parameters.time_diff_matrix,
            model_parameters.subject_birth_ix[ix_subject],
            AbstractArray{Bool}(view(infections, :, ix_subject)),
            obs_lookup[ix_subject], 
            y_pred
        )

         observed_titre[ix_subject] ~ TitreArrayNormal(y_pred, obs_sd, obs_min, obs_max)
    end
end

function create_special_matrix_iterative(n::Int)
    if n <= 0
        throw(ArgumentError("Matrix size n must be a positive integer."))
    end

    M = zeros(Float64, n, n)
    for i in 1:n
        for j in 1:n
            if i == j       # Main diagonal
                M[i, j] = 1.0
            elseif abs(i - j) == 1 # Super-diagonal (j = i+1) or Sub-diagonal (j = i-1)
                M[i, j] = 0.5
            end
            # Other elements remain 0.0 by initialization
        end
    end
    return M
end


function ppd_infer_distances(chain, model_parameters, n_ppd_subjects, n_draws)
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
    ix_strain_gaps = findall(s -> startswith(s, "strain_gaps"), col_names)

    for ix_draw in 1:n_draws
        ix_sample = sample(1:nrow(chain)) 

        draw = chain[ix_sample,:]
        
        # TODO cache across ix_draw
        obs_lookup = make_obs_lookup(obs_df_grouped[ix_draw])
        obs_view = make_obs_views(obs_df_grouped[ix_draw])

        infections = reshape(Vector(chain[ix_sample, ix_infections]), model_parameters.n_t_steps, model_parameters.n_subjects)

        strain_gaps = Vector(chain[ix_sample, ix_strain_gaps])
        strain_locations = vcat([0], cumsum(exp.(strain_gaps)))
    
        antigenic_distances = zeros(model_parameters.n_t_steps, model_parameters.n_t_steps)

        for j in 1:model_parameters.n_t_steps
            for i in 1:j
                loc_i = strain_locations[i]
                loc_j = strain_locations[j]

                antigenic_distances[i, j] = sqrt((loc_i - loc_j)^2)
                antigenic_distances[j, i] = antigenic_distances[i, j]
            end
        end

        for ix_subject in 1:n_ppd_subjects
            waning_curve_individual!(
                draw.mu_long, draw.mu_short, draw.omega,
                0.13, 0.03, draw.tau,

                antigenic_distances, model_parameters.time_diff_matrix,
                model_parameters.subject_birth_ix[ix_subject],

                AbstractArray{Bool}(view(infections, :, ix_subject)),

                obs_lookup[ix_subject], 
                view(obs_df_grouped[ix_draw].observed_titre, obs_view[ix_subject])
            )
        end
    end

    return DataFrame(obs_df_grouped)
end