

struct GibbsExplorer3
    proposal_function::Function
    symbols_not_inf::Vector{Symbol}
    n_repeats::Int

    sp::StaticModelParameters

    sigma_covar::Float64
end

function Pigeons.step!(explorer::GibbsExplorer3, replica, shared)
    state = replica.state
    rng = replica.rng
    recorders = replica.recorders
    chain = replica.chain

    log_potential = Pigeons.find_log_potential(replica, shared.tempering, shared)

    step_slicing(state, explorer, log_potential, replica)
    for i in explorer.n_repeats
        # step_rmh(state, explorer, rng, log_potential, recorders, chain)
        step_infections(state, explorer, rng, log_potential)
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

function step_rmh(state, explorer, rng, log_potential, recorders, chain)
    # Add zero to initialise
    Pigeons.@record_if_requested!(recorders, :n_accepted, (chain, 0))
    Pigeons.@record_if_requested!(recorders, :n_rejected, (chain, 0)) 

    # Get the logdensity function
    logprob_previous = log_potential(state)

    theta_previous = zeros(length(explorer.symbols_not_inf))
    for (i, meta) in enumerate(state.metadata[explorer.symbols_not_inf])
        theta_previous[i] = meta.vals[1]
    end


    # Theta proposal
    # Note no log-transform here as bijection is provided by Pigeons.jl
    theta_new_dist = MvNormal(theta_previous, explorer.sigma_covar)
    theta_new = rand(rng, theta_new_dist)

    for (i, meta) in enumerate(state.metadata[explorer.symbols_not_inf])
        meta.vals[1] = theta_new[i]
    end

    logprob_proposal = log_potential(state)

    log_target_ratio = logprob_proposal - logprob_previous

    if -Random.randexp(rng) <= log_target_ratio
        Pigeons.@record_if_requested!(recorders, :n_accepted, (chain, 1)) 
    else   
        Pigeons.@record_if_requested!(recorders, :n_rejected, (chain, 1)) 
        for (i, meta) in enumerate(state.metadata[explorer.symbols_not_inf])
            meta.vals[1] = theta_previous[i]
        end
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


function Pigeons.adapt_explorer(
    explorer::GibbsExplorer3,
    reduced_recorders,
    current_pt,
    new_tempering
)
    # n_accepted = Pigeons.value.(values(Pigeons.value(reduced_recorders.n_accepted)))
    # n_rejected = Pigeons.value.(values(Pigeons.value(reduced_recorders.n_rejected)))

    # n_steps = n_accepted[1] + n_rejected[1] * explorer.n_repeats

    # pr_accept = mean(n_accepted ./ (n_accepted .+ n_rejected))

    # This is a bit weird
    # because this adaptation only happens once (rather than each iteration)
    # but does seem to work OK?
    # sigma_covar_adapted = exp(log(explorer.sigma_covar) + (pr_accept - 0.234) * 0.999 ^ n_steps)

    # @printf("pr_accept = %6.4f; sigma_covar_adapted = %6.4f;\n", pr_accept, sigma_covar_adapted)

    sigma_covar_adapted = 0.01

    return GibbsExplorer3(
        explorer.proposal_function,
        explorer.symbols_not_inf,
        explorer.n_repeats,
        explorer.sp,
        sigma_covar_adapted
    )
end

n_accepted() = Pigeons.GroupBy(Int, Pigeons.Sum())
n_rejected() = Pigeons.GroupBy(Int, Pigeons.Sum())

function Pigeons.explorer_recorder_builders(explorer::GibbsExplorer3)
    result = [
        n_accepted,
        n_rejected
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