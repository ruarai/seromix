
# The Kucharski (2018) model for waning immunity
@model function waning_model_kucharski(
    model_parameters::FixedModelParameters,
    prior_infection_dist::Distribution,
    use_corrected_titre::Bool,

    obs_lookup_strain, obs_lookup_ix, obs_views,
    n_max_ind_obs::Int,

    observed_titre::Vector{Vector{Float64}};

    mixture_importance_sampling = false
)
    mu_long ~ Uniform(0.0, 10.0)
    mu_short ~ Uniform(0.0, 10.0)


    omega ~ Uniform(0.0, 1.0)

    sigma_long ~ Uniform(0.0, 10.0)
    sigma_short ~ Uniform(0.0, 10.0)

    tau ~ Uniform(0.0, 10.0)

    obs_sd ~ Uniform(0.0, 10.0)

    infections ~ prior_infection_dist

    # If we're in an "IndividualSubsetContext", only calculate the likelihood over a single subject
    # Otherwise, calculate across all subjects
    subjects_to_process = if DynamicPPL.leafcontext(__context__) isa IndividualSubsetContext
        SA[DynamicPPL.leafcontext(__context__).ix] # StaticArray with one element (context.ix)
    else
        1:model_parameters.n_subjects
    end

    # Reduce the total memory allocation across the observed titre
    # by calculating across a single pre-allocated array.
    y_pred_buffer = Vector{Float64}(undef, n_max_ind_obs)

    # We calculate streaming logsumexp(-logp) to add to the likelihood
    # at the end if doing mixture importance sampling (otherwise, we add zero)
    # https://www.nowozin.net/sebastian/blog/streaming-log-sum-exp-computation.html
    r, alpha = mixture_importance_sampling ? (0.0, -Inf) : (1.0, 0.0)

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

        if mixture_importance_sampling
            lp = 0.0
            # TODO move to function which modifies r, alpha, somehow
            @inbounds for i in 1:length(y_pred)
                lpp = titre_logpdf_component(
                    observed_titre[ix_subject][i],
                    y_pred[i],
                    obs_sd, const_titre_min, const_titre_max
                )

                lp += lpp

                neg_lp = -lpp
                r = neg_lp <= alpha ? r + exp(neg_lp - alpha) : r * exp(alpha - neg_lp) + 1.0
                alpha = neg_lp <= alpha ? alpha : neg_lp
            end

            @addlogprob! lp
        else
            observed_titre[ix_subject] ~ TitreArrayNormal(y_pred, obs_sd, const_titre_min, const_titre_max, use_corrected_titre)
        end

    end

    log_sum_exp_neg_lp = log(r) + alpha
    @addlogprob! log_sum_exp_neg_lp
end


