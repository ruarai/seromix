
struct FixedModelParameters
    n_t_steps::Int
    n_subjects::Int

    antigenic_distances::Matrix{Float64}

    time_diff_matrix::Matrix{Float64}

    # The first time-step in which the subject may be infected
    # if zero, assume that they were born prior to simulation,
    # such that it is equivalent to a value of one
    subject_birth_ix::Vector{Int}
end

struct WaningModelCache
    obs_lookup_strain::Vector{Dict{Int64, Vector{Int64}}}
    obs_lookup_ix::Vector{Dict{Int64, Vector{Int64}}}

    obs_views::Vector{UnitRange{Int64}}

    n_max_ind_obs::Int
end

function make_waning_model(
    model_parameters::FixedModelParameters,
    obs_df::DataFrame;
    prior_infection_dist::Distribution,
    turing_model = waning_model_kucharski,
    mixture_importance_sampling = false,
    use_corrected_titre = true
)
    all(diff(obs_df.ix_subject) .>= 0) || throw(ArgumentError("ix_subject in obs_df must be sorted in ascending order."))

    obs_views = make_obs_views(obs_df)

    individual_titre_obs = [obs_df.observed_titre[v] for v in obs_views]

    obs_lookup_strain, obs_lookup_ix = make_obs_lookup(obs_df)
    n_max_ind_obs = maximum(length.(obs_views))

    model_cache = WaningModelCache(
        obs_lookup_strain, obs_lookup_ix,
        obs_views,
        n_max_ind_obs
    )
    
    return turing_model(
        model_parameters,
        prior_infection_dist,
        individual_titre_obs,
        model_cache;
        
        mixture_importance_sampling = mixture_importance_sampling,
        use_corrected_titre = use_corrected_titre
    );
end


function general_waning_likelihood(
    params, infections, model_parameters,
    observed_titre::Vector{Vector{Float64}},
    model_cache::WaningModelCache,
    individual_waning_function::Function,
    context::DynamicPPL.AbstractContext;
    mixture_importance_sampling::Bool = false,
    use_corrected_titre::Bool = true
)
    # If we're in an "IndividualSubsetContext", only calculate the likelihood over a single subject
    # Otherwise, calculate across all subjects
    subjects_to_process = if context isa IndividualSubsetContext
        SA[context.ix] # StaticArray with one element (context.ix)
    else
        1:model_parameters.n_subjects
    end

    # Reduce the total memory allocation across the observed titre
    # by calculating across a single pre-allocated array.
    y_pred_buffer = Vector{Float64}(undef, model_cache.n_max_ind_obs)

    # We calculate streaming logsumexp(-logp) to add to the likelihood
    # at the end if doing mixture importance sampling (otherwise, we add zero), see:
    # https://mc-stan.org/loo/articles/loo2-mixis.html
    # https://www.nowozin.net/sebastian/blog/streaming-log-sum-exp-computation.html
    r, alpha = (0.0, -Inf)

    # Accumate the log-likelihood
    lp = 0.0

    for ix_subject in subjects_to_process
        n_obs_subject = length(model_cache.obs_views[ix_subject])
        y_pred = view(y_pred_buffer, 1:n_obs_subject)
        fill!(y_pred, 0.0)

        # Fills in y_pred according to antibody dynamics
        individual_waning_function(
            params,
            model_parameters.antigenic_distances, model_parameters.time_diff_matrix,
            model_parameters.subject_birth_ix[ix_subject],
            view(infections, :, ix_subject),
            model_cache.obs_lookup_strain[ix_subject],
            model_cache.obs_lookup_ix[ix_subject],
            y_pred
        )

        # Calculate the log-likelihood of y_pred for this ix_subject
        if mixture_importance_sampling
            @inbounds for i in 1:length(y_pred)
                # Pointwise log likelihood
                lpp = if use_corrected_titre
                    titre_logpdf_component(
                        observed_titre[ix_subject][i],
                        y_pred[i],
                        params.obs_sd, const_titre_min, const_titre_max
                    )
                else
                    titre_logpdf_component_uncorrected(
                        observed_titre[ix_subject][i],
                        y_pred[i],
                        params.obs_sd, const_titre_min, const_titre_max
                    )
                end

                lp += lpp

                # Accumulation for the logsumexp(-lpp[ix_subject])
                neg_lpp = -lpp
                r = neg_lpp <= alpha ? r + exp(neg_lpp - alpha) : r * exp(alpha - neg_lpp) + 1.0
                alpha = neg_lpp <= alpha ? alpha : neg_lpp
            end
        else
            if use_corrected_titre
                lp += apply_logpdf_simd(observed_titre[ix_subject], y_pred, params.obs_sd, const_titre_min, const_titre_max)
            else
                lp += apply_logpdf_uncorrected(observed_titre[ix_subject], y_pred, params.obs_sd, const_titre_min, const_titre_max)
            end
        end

    end

    if mixture_importance_sampling
        lp += log(r) + alpha
    end

    return lp
end

function general_waning_likelihood(
    params, infections::Matrix{Float64}, model_parameters,
    observed_titre::Vector{Vector{Float64}},
    model_cache::WaningModelCache,
    individual_waning_function::Function,
    context::DynamicPPL.AbstractContext;
    mixture_importance_sampling::Bool = false,
    use_corrected_titre::Bool = true
)
    return general_waning_likelihood(
        params, convert.(Bool, infections), model_parameters,
        observed_titre,
        model_cache,
        individual_waning_function,
        context;
        mixture_importance_sampling = mixture_importance_sampling,
        use_corrected_titre = use_corrected_titre
    )
end

function waning_curve!(
    params, individual_waning_function,
    dist_matrix::Matrix{Float64},
    time_diff_matrix::Matrix{Float64},
    subject_birth_ix::Vector{Int},
    infections::Matrix{Bool},
    obs_lookup_strain,
    obs_lookup_ix,
    obs_views::Vector{UnitRange{Int}},
    y::AbstractArray{T}
) where T <: Real
    for ix_subject in axes(infections, 2)
        individual_waning_function(
            params,

            dist_matrix, time_diff_matrix,
            subject_birth_ix[ix_subject],

            view(infections, :, ix_subject),
            obs_lookup_strain[ix_subject],
            obs_lookup_ix[ix_subject],

            view(y, obs_views[ix_subject])
        )
    end
    
    return y
end



function pointwise_likelihood(
    chain_df, model_data,
    turing_model
)
    p = read_model_parameters(model_data)

    obs_df = DataFrame(model_data["observations"])

    col_names = names(chain_df)
    ix_infections = findall(s -> startswith(s, "infections"), col_names)

    obs_lookup_strain, obs_lookup_ix = make_obs_lookup(obs_df)
    obs_views = make_obs_views(obs_df)

    logp = zeros(nrow(chain_df), nrow(obs_df))

    param_symbols = model_symbols_apart_from(turing_model, [:infections])

    individual_waning_function = turing_function_to_waning_function(model.f)


    y_pred_buffer = zeros(maximum([length(v) for v in obs_views]))

    @showprogress for ix_sample in 1:nrow(chain_df)
        draw = chain_df[ix_sample, :]
        

        infections = reshape(Vector(chain_df[ix_sample, ix_infections]), p.n_t_steps, p.n_subjects)

        for ix_subject in 1:p.n_subjects
            n_obs_subject = length(obs_views[ix_subject])
            y_pred = view(y_pred_buffer, 1:n_obs_subject)
            fill!(y_pred, 0.0)

            params = NamedTuple{param_symbols}([draw[i] for i in param_symbols])

            individual_waning_function(
                params,

                p.antigenic_distances, p.time_diff_matrix,
                p.subject_birth_ix[ix_subject],

                AbstractArray{Bool}(view(infections, :, ix_subject)),

                obs_lookup_strain[ix_subject], obs_lookup_ix[ix_subject], 
                y_pred
            )

            @inbounds for ix_obs in 1:n_obs_subject
                ix_obs_absolute = obs_views[ix_subject][ix_obs]
                logp[ix_sample, ix_obs_absolute] = titre_logpdf_component(
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