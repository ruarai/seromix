


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
    param_link_dists = [
        Uniform(0, 10), Uniform(0, 10),
        Uniform(0, 1), Uniform(0, 10),
        Uniform(0, 10), Uniform(0, 10),
        Uniform(0, 10)
    ]
    param_means = [2.0, 2.5, 0.8, 0.15, 0.05, 0.05, 1.0]
    means = [link(param_link_dists[i], param_means[i]) for i in 1:7]

    params ~ MvNormal(means, I)

    mu_long := invlink(param_link_dists[1], params[1])
    mu_short := invlink(param_link_dists[2], params[2])
    omega := invlink(param_link_dists[3], params[3])
    sigma_long := invlink(param_link_dists[4], params[4])
    sigma_short := invlink(param_link_dists[5], params[5])
    tau := invlink(param_link_dists[6], params[6])
    obs_sd := invlink(param_link_dists[7], params[7])
    
    obs_min = convert(typeof(mu_long), const_titre_min)
    obs_max = convert(typeof(mu_long), const_titre_max)

    row_means ~ MvNormal(fill(-1.0, p.n_t_steps), I)

    infections ~ MatrixHierarchicalBernoulli(
        row_means,
        zeros(typeof(mu_long), p.n_subjects),
        p.n_t_steps,
        p.n_subjects
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