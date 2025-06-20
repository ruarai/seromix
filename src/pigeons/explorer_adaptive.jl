

struct StateExplorerAdaptive2
    # slice_sampler
    proposal_function::Function

    symbols_not_inf::Vector{Symbol}

    n_t_steps::Int
    n_subjects::Int

    sigma_step::Float64
    inf_step::Float64
end


function Pigeons.step!(explorer::StateExplorerAdaptive2, replica, shared)
    state = replica.state
    rng = replica.rng
    recorders = replica.recorders
    chain = replica.chain

    log_potential = Pigeons.find_log_potential(replica, shared.tempering, shared)

    # Slice sampling
    h = SliceSampler()
    cached_lp = -Inf
    for meta in state.metadata[explorer.symbols_not_inf]
        cached_lp = Pigeons.slice_sample!(h, meta.vals, log_potential, cached_lp, replica)
    end


    # Infections

    n_accept = 0
    n_reject = 0

    inf_vi = state.metadata.infections
    
    n_t_steps = explorer.n_t_steps
    n_subjects = explorer.n_subjects

    for ix_subject in 1:n_subjects
        swap_indices, log_hastings_ratio = explorer.proposal_function(
            rng,
            inf_vi.vals,
            ix_subject,
            n_t_steps,
            explorer.inf_step
        )

        log_pr_before = log_potential(state, IndividualSubsetContext(ix_subject))

        apply_swaps!(inf_vi.vals, swap_indices, ix_subject, n_t_steps)

        log_pr_after = log_potential(state, IndividualSubsetContext(ix_subject))

        accept_ratio = exp(log_pr_after - log_pr_before + log_hastings_ratio)
        if accept_ratio < 1 && rand(rng) > accept_ratio
            # reject: revert the move we just proposed
            apply_swaps!(inf_vi.vals, swap_indices, ix_subject, n_t_steps)

            n_reject += 1
        else
            n_accept += 1
        end
    end

    Pigeons.@record_if_requested!(recorders, :n_accepted_inf, (chain, n_accept)) 
    Pigeons.@record_if_requested!(recorders, :n_rejected_inf, (chain, n_reject)) 
end

function Pigeons.adapt_explorer(
    explorer::StateExplorerAdaptive2,
    reduced_recorders,
    current_pt,
    new_tempering
)
    # Adjust random walk
    mean_accepted_param = mean(Pigeons.value.(values(Pigeons.value(reduced_recorders.n_accepted_param))))
    mean_rejected_param = mean(Pigeons.value.(values(Pigeons.value(reduced_recorders.n_rejected_param))))

    target_accept_rate = 0.234

    pr_accept_param = mean_accepted_param / (mean_accepted_param + mean_rejected_param)

    sigma_step_adapted = clamp(exp(log(explorer.sigma_step) + pr_accept_param - target_accept_rate), 1e-5, 1.0)

    println("Param acceptance rate = $(round(pr_accept_param, digits = 2)), adapted step size = $(round(sigma_step_adapted, digits = 2))")

    mean_accepted_inf = mean(Pigeons.value.(values(Pigeons.value(reduced_recorders.n_accepted_inf))))
    mean_rejected_inf = mean(Pigeons.value.(values(Pigeons.value(reduced_recorders.n_rejected_inf))))

    pr_accept_inf = mean_accepted_inf / (mean_accepted_inf + mean_rejected_inf)

    inf_step_adapted = clamp(exp(log(explorer.inf_step) + pr_accept_inf - target_accept_rate), 1e-5, 10.0)

    println("Inf acceptance rate = $(round(pr_accept_inf, digits = 2)), adapted inf step size = $(round(inf_step_adapted, digits = 2))")


    return StateExplorerAdaptive2(
        explorer.proposal_function,
        explorer.symbols_not_inf,
        explorer.n_t_steps,
        explorer.n_subjects,
        sigma_step_adapted,
        inf_step_adapted
    )
end


n_accepted_param() = Pigeons.GroupBy(Int, Pigeons.Sum())
n_rejected_param() = Pigeons.GroupBy(Int, Pigeons.Sum())
n_accepted_inf() = Pigeons.GroupBy(Int, Pigeons.Sum())
n_rejected_inf() = Pigeons.GroupBy(Int, Pigeons.Sum())

function Pigeons.explorer_recorder_builders(explorer::StateExplorerAdaptive2)
    result = [
        n_accepted_param,
        n_rejected_param,
        n_accepted_inf,
        n_rejected_inf
    ]
    return result
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