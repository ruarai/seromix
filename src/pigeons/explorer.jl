

struct GibbsExplorer4
    proposal_function::Function
    symbols_not_inf::Vector{Symbol}
    n_repeats::Int

    sp::StaticModelParameters
end

function Pigeons.step!(explorer::GibbsExplorer4, replica, shared)
    state = replica.state
    rng = replica.rng
    recorders = replica.recorders
    chain = replica.chain

    log_potential = Pigeons.find_log_potential(replica, shared.tempering, shared)


    for _ in 1:explorer.n_repeats
        step_slicing(state, explorer, log_potential, replica)
        for _ in 1:explorer.n_repeats
            step_infections(state, explorer, rng, log_potential)
        end
    end
end

function step_slicing(state, explorer, log_potential, replica)
    h = SliceSampler()

    # Slice sampling
    cached_lp = -Inf
    for meta in state.metadata[explorer.symbols_not_inf]
        cached_lp = Pigeons.slice_sample!(h, meta.vals, log_potential, cached_lp, replica)
    end
end

function step_infections(state, explorer, rng, log_potential)
    inf_vi = state.metadata.infections
    
    n_t_steps = explorer.sp.n_t_steps
    n_subjects = explorer.sp.n_subjects

    subject_indices = sample(rng, 1:n_subjects, ceil(Int, 0.4 * n_subjects), replace = false)
    for ix_subject in subject_indices
        swap_indices, log_hastings_ratio = explorer.proposal_function(
            rng,
            inf_vi.vals,
            ix_subject,
            n_t_steps,
            explorer.sp.subject_birth_ix[ix_subject]
        )

        log_pr_before = log_potential(state, IndividualSubsetContext(ix_subject))

        apply_swaps!(inf_vi.vals, swap_indices, ix_subject, n_t_steps, explorer.sp.subject_birth_ix[ix_subject])

        log_pr_after = log_potential(state, IndividualSubsetContext(ix_subject))

        accept_ratio = exp(log_pr_after - log_pr_before + log_hastings_ratio)
        if accept_ratio < 1 && rand(rng) > accept_ratio
            # reject: revert the move we just proposed
            apply_swaps!(inf_vi.vals, swap_indices, ix_subject, n_t_steps, explorer.sp.subject_birth_ix[ix_subject])
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


Bijectors.bijector(d::AbstractMatrixBernoulli) = identity