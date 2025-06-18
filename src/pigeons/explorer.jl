

struct CustomExplorer2
    slice_sampler

    symbols_not_inf

    n_t_steps
    n_subjects
end


function Pigeons.step!(explorer::CustomExplorer2, replica, shared)
    state = replica.state
    rng = replica.rng
    h = explorer.slice_sampler

    log_potential = Pigeons.find_log_potential(replica, shared.tempering, shared)

    # log_theta_current = zeros(length(explorer.symbols_not_inf))

    # for i in eachindex(explorer.symbols_not_inf)
    #     meta = state.metadata[explorer.symbols_not_inf[i]]
    #     log_theta_current[i] = meta.vals[1]
    # end

    # logprob_previous = log_potential(state)

    # log_theta_new_dist = MvNormal(log_theta_current, I * 1.0)
    # log_theta_new = rand(rng, log_theta_new_dist)

    # for i in eachindex(explorer.symbols_not_inf)
    #     meta = state.metadata[explorer.symbols_not_inf[i]]
    #     meta.vals[1] = log_theta_new[i]
    # end

    # logprob_proposal = log_potential(state)

    # log_target_ratio = logprob_proposal - logprob_previous

    # if Random.randexp(rng) > log_target_ratio
    #     # Reject, revert changes
    #     for i in eachindex(explorer.symbols_not_inf)
    #         meta = state.metadata[explorer.symbols_not_inf[i]]
    #         meta.vals[1] = log_theta_current[i]
    #     end
    # end

    cached_lp = -Inf
    for meta in state.metadata[explorer.symbols_not_inf]
        cached_lp = Pigeons.slice_sample!(h, meta.vals, log_potential, cached_lp, replica)
    end


    # Infections
    inf_vi = state.metadata.infections
    
    n_t_steps = explorer.n_t_steps
    n_subjects = explorer.n_subjects

    for ix_subject in 1:n_subjects
        swap_indices, log_hastings_ratio = proposal_original_corrected(
            rng,
            inf_vi.vals,
            ix_subject,
            n_t_steps
        )

        log_pr_before = log_potential(state, IndividualSubsetContext(ix_subject))

        apply_swaps!(inf_vi.vals, swap_indices, ix_subject, n_t_steps)

        log_pr_after = log_potential(state, IndividualSubsetContext(ix_subject))

        accept_ratio = exp(log_pr_after - log_pr_before + log_hastings_ratio)
        if accept_ratio < 1 && rand(rng) > accept_ratio
            # reject: revert the move we just proposed
            apply_swaps!(inf_vi.vals, swap_indices, ix_subject, n_t_steps)
        end
    end
end


function chain_infections_matrix_2(chain, ix_iter, ix_chain, params)
    infections = zeros(Bool, params.n_t_steps, params.n_subjects)

    for ix_subject in 1:params.n_subjects
        for ix_t in max(1, params.subject_birth_ix[ix_subject]):params.n_t_steps
            ix = (ix_subject - 1) * params.n_t_steps + ix_t
            infections[ix_t, ix_subject] = chain[ix_iter, Symbol("infections[$ix]"), ix_chain]
        end
    end

    return infections
end


function chain_infections_prob_2(chain, params)
    infections = zeros(Float64, params.n_t_steps, params.n_subjects)

    for ix_subject in 1:params.n_subjects
        for ix_t in max(1, params.subject_birth_ix[ix_subject]):params.n_t_steps
            ix = (ix_subject - 1) * params.n_t_steps + ix_t
            infections[ix_t, ix_subject] = mean(chain[Symbol("infections[$ix]")])
        end
    end

    return infections
end


(log_potential::Pigeons.TuringLogPotential{<:Any,<:DynamicPPL.DefaultContext})(vi, context) =
    try

        return DynamicPPL.getlogp(last(DynamicPPL.evaluate!!(log_potential.model, vi, context)))
    catch e
        (isa(e, DomainError) || isa(e, BoundsError)) && return -Inf
        rethrow(e)
    end


(interpolated::Pigeons.InterpolatedLogPotential)(x, context) = 
    if interpolated.beta == zero(interpolated.beta)
        interpolated.path.ref(x)
    elseif interpolated.beta == one(interpolated.beta)
        interpolated.path.target(x, context)
    else
        Pigeons.interpolate(interpolated.path.interpolator, interpolated.path.ref(x), interpolated.path.target(x, context), interpolated.beta)
    end