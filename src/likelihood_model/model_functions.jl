
struct StaticModelParameters
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

function WaningModelCache(obs_df)
    obs_lookup_strain, obs_lookup_ix = make_obs_lookup(obs_df)
    n_max_ind_obs = maximum(length.(obs_views))

    model_cache = WaningModelCache(
        obs_lookup_strain, obs_lookup_ix,
        obs_views,
        n_max_ind_obs
    )
    
    return model_cache
end

function make_waning_model(
    sp::StaticModelParameters,
    obs_df::DataFrame;
    prior_infection_dist::Distribution,
    turing_model = waning_model_kucharski,
    mixture_importance_sampling = false,
    use_corrected_titre = true
)
    all(diff(obs_df.ix_subject) .>= 0) || throw(ArgumentError("ix_subject in obs_df must be sorted in ascending order."))

    obs_views = make_obs_views(obs_df)

    individual_titre_obs = [obs_df.observed_titre[v] for v in obs_views]

    model_cache = WaningModelCache(obs_df)
    
    return turing_model(
        sp,
        prior_infection_dist,
        individual_titre_obs,
        model_cache;
        
        mixture_importance_sampling = mixture_importance_sampling,
        use_corrected_titre = use_corrected_titre
    );
end

function general_waning_likelihood(
    params,
    infections,
    observed_titre::Vector{Vector{Float64}},

    sp::StaticModelParameters,
    model_cache::WaningModelCache,
    __context__,

    individual_waning_function::Function;

    use_corrected_titre::Bool = true,
    mixture_importance_sampling::Bool = false
)
    leaf_context = DynamicPPL.leafcontext(__context__)
    if leaf_context isa DynamicPPL.PriorContext
        return 0.0
    end

    # If we're in an "IndividualSubsetContext", only calculate the likelihood over a single subject
    # Otherwise, calculate across all subjects
    subjects_to_process = if leaf_context isa IndividualSubsetContext
        SA[leaf_context.ix] # StaticArray with one element (leaf_context.ix)
    else
        1:sp.n_subjects
    end

    # Reduce the total memory allocation across the observed titre
    # by calculating across a single pre-allocated array.
    latent_titre_buffer = Vector{typeof(params[1])}(undef, model_cache.n_max_ind_obs)

    # We calculate streaming logsumexp(-logp) to add to the likelihood
    # at the end if doing mixture importance sampling (otherwise, we add zero), see:
    # https://mc-stan.org/loo/articles/loo2-mixis.html
    # https://www.nowozin.net/sebastian/blog/streaming-log-sum-exp-computation.html
    r, alpha = (0.0, -Inf)

    # Accumate the log-likelihood
    lp = 0.0

    for ix_subject in subjects_to_process
        n_obs_subject = length(model_cache.obs_views[ix_subject])
        latent_titre = view(latent_titre_buffer, 1:n_obs_subject)
        fill!(latent_titre, 0.0)

        # Fills in y_pred according to antibody dynamics
        individual_waning_function(
            params,
            view(infections, :, ix_subject),
            latent_titre,
            ix_subject,

            sp,
            model_cache
        )

        # Calculate the log-likelihood of y_pred for this ix_subject
        if mixture_importance_sampling
            @inbounds for i in 1:length(latent_titre)
                # Pointwise log likelihood
                lpp = if use_corrected_titre
                    titre_logpdf_component(
                        observed_titre[ix_subject][i],
                        latent_titre[i],
                        params.obs_sd, const_titre_min, const_titre_max
                    )
                else
                    titre_logpdf_component_uncorrected(
                        observed_titre[ix_subject][i],
                        latent_titre[i],
                        params.obs_sd, const_titre_min, const_titre_max
                    )
                end

                if isinf(lpp) && lpp < 0
                    @show lpp
                    @show i, ix_subject
                    @show observed_titre[ix_subject][i]
                    @show latent_titre[i]
                    @show params.obs_sd
                    # It might be useful to see the state of r and alpha too
                    @show r, alpha
                end

                lp += lpp

                # Accumulation for the logsumexp(-lpp[ix_subject])
                neg_lpp = -lpp
                r = neg_lpp <= alpha ? r + exp(neg_lpp - alpha) : r * exp(alpha - neg_lpp) + 1.0
                alpha = neg_lpp <= alpha ? alpha : neg_lpp
            end
        else
            if use_corrected_titre
                lp += apply_logpdf_simd(observed_titre[ix_subject], latent_titre, params.obs_sd, const_titre_min, const_titre_max)
            else
                lp += apply_logpdf_uncorrected(observed_titre[ix_subject], latent_titre, params.obs_sd, const_titre_min, const_titre_max)
            end
        end

    end

    if mixture_importance_sampling
        lp += log(r) + alpha
    end

    return lp
end

function general_waning_likelihood(
    params,
    infections::Matrix{Float64},
    observed_titre::Vector{Vector{Float64}},

    sp::StaticModelParameters,
    model_cache::WaningModelCache,
    __context__,

    individual_waning_function::Function;

    use_corrected_titre::Bool = true,
    mixture_importance_sampling::Bool = false
)
    return general_waning_likelihood(
        params, convert.(Bool, infections), observed_titre, sp, model_cache, __context__, individual_waning_function;
        use_corrected_titre = use_corrected_titre, mixture_importance_sampling = mixture_importance_sampling
    )
end

function waning_curve!(
    params,
    individual_waning_function,
    sp::StaticModelParameters,
    infections::Matrix{Bool},
    model_cache::WaningModelCache,
    latent_titre::AbstractArray{T}
) where T <: Real
    for ix_subject in axes(infections, 2)
        individual_waning_function(
            params,
            view(infections, :, ix_subject),
            latent_titre,
            ix_subject,

            sp,
            model_cache
        )
    end
end



function model_pointwise_likelihood(
    chain_df, sp, obs_df,
    turing_model;
    fixed_params = missing
)
    col_names = names(chain_df)
    ix_infections = findall(s -> startswith(s, "infections"), col_names)

    model_cache = WaningModelCache(obs_df)

    logp = zeros(nrow(chain_df), nrow(obs_df))

    param_symbols = model_symbols_apart_from(turing_model, [:infections])

    individual_waning_function = turing_function_to_waning_function(turing_model.f)

    @showprogress Threads.@threads for ix_sample in 1:nrow(chain_df)
        draw = chain_df[ix_sample, :]
        latent_titre_buffer = zeros(model_cache.n_max_ind_obs)
        

        infections = reshape(Vector(chain_df[ix_sample, ix_infections]),sp.n_t_steps,sp.n_subjects)
        infections = convert(Matrix{Bool}, convert.(Bool, infections))

        
        params = NamedTuple{param_symbols}([draw[i] for i in param_symbols])
        
        if !ismissing(fixed_params)
            params = merge(params, fixed_params)
        end

        for ix_subject in 1:sp.n_subjects
            n_obs_subject = length(obs_views[ix_subject])
            latent_titre = view(latent_titre_buffer, 1:n_obs_subject)
            fill!(latent_titre, 0.0)

            individual_waning_function(
                params,
                AbstractArray{Bool}(view(infections, :, ix_subject)),
                latent_titre,
                ix_subject,

                sp,
                model_cache
            )

            @inbounds for ix_obs in 1:n_obs_subject
                ix_obs_absolute = obs_views[ix_subject][ix_obs]
                logp[ix_sample, ix_obs_absolute] = titre_logpdf_component(
                    obs_df.observed_titre[ix_obs_absolute],
                    latent_titre[ix_obs],
                    draw.obs_sd,
                    const_titre_min,
                    const_titre_max
                )
            end
        end
    end

    return logp
end

function model_sum_mixIS(chain_df, sp, obs_df, model)
    chains = unique(chain_df.chain)
    elpd_mixIS_sum = zeros(length(chains))

    lpp = model_pointwise_likelihood(chain_df, sp, obs_df, model)

    for ix_chain in chains
        rows_chain = chain_df.chain .== ix_chain

        l_common_mix = logsumexp(-lpp[rows_chain,:], dims = 2)
        log_weights = -lpp[rows_chain,:] .- l_common_mix
        elpd_mixIS = logsumexp(-l_common_mix) .- logsumexp(log_weights, dims = 1)
        elpd_mixIS_sum[ix_chain] = sum(elpd_mixIS)
    end

    return elpd_mixIS_sum
end


function model_ppd(
    chain_df, sp,
    turing_model;
    n_samples = 10,
    fixed_params = missing
)
    samples_indices = sort(sample(1:size(chain_df, 1), n_samples, replace = false))

    complete_obs = expand_grid(
        ix_t_obs = 1:sp.n_t_steps,
        ix_strain = 1:sp.n_t_steps,
        ix_subject = 1:sp.n_subjects,
        observed_titre = 0.0,
        ix_sample = eachindex(samples_indices)
    )

    complete_obs_grouped = groupby(complete_obs, :ix_sample)


    col_names = names(chain_df)
    ix_infections = findall(s -> startswith(s, "infections"), col_names)


    param_symbols = model_symbols_apart_from(turing_model, [:infections])

    individual_waning_function = turing_function_to_waning_function(turing_model.f)

    @showprogress Threads.@threads for ix_sample in eachindex(samples_indices)

        obs_sample = complete_obs_grouped[(ix_sample = ix_sample, )]


        model_cache = WaningModelCache(obs_sample)


        draw = chain_df[samples_indices[ix_sample], :]


        infections = reshape(Vector(chain_df[samples_indices[ix_sample], ix_infections]),sp.n_t_steps,sp.n_subjects)
        infections = convert(Matrix{Bool}, convert.(Bool, infections))


        params = NamedTuple{param_symbols}([draw[i] for i in param_symbols])

        if !ismissing(fixed_params)
            params = merge(params, fixed_params)
        end

        obs_val = zeros(nrow(obs_sample))
        
        waning_curve!(
            params,
            individual_waning_function,
            sp,
            infections,
            model_cache,
            obs_val
        )

        complete_obs.observed_titre[complete_obs.ix_sample .== ix_sample] .= obs_val
    end

    return complete_obs
end