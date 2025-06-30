
using AdvancedMH

# Taken from RobustAdaptiveMetropolis in
# https://github.com/TuringLang/AdvancedMH.jl
# adapted to disable the warmup steps
# (i.e. is always adapting)

Base.@kwdef struct RobustAdaptiveMetropolis2{T,A<:Union{Nothing,AbstractMatrix{T}}} <:
                   AdvancedMH.MHSampler
    "target acceptance rate. Default: 0.234."
    α::T = 0.234
    "negative exponent of the adaptation decay rate. Default: `0.6`."
    γ::T = 0.6
    "initial lower-triangular Cholesky factor of the covariance matrix. If specified, should be convertible into a `LowerTriangular`. Default: `nothing`, which is interpreted as the identity matrix."
    S::A = nothing
    "lower bound on eigenvalues of the adapted Cholesky factor. Default: `0.0`."
    eigenvalue_lower_bound::T = 0.0
    "upper bound on eigenvalues of the adapted Cholesky factor. Default: `Inf`."
    eigenvalue_upper_bound::T = Inf
end

struct RobustAdaptiveMetropolis2State{T1,L,A,T2,T3}
    "current realization of the chain."
    x::T1
    "log density of `x` under the target model."
    logprob::L
    "current lower-triangular Cholesky factor."
    S::A
    "log acceptance ratio of the previous iteration (not necessarily of `x`)."
    logα::T2
    "current step size for adaptation of `S`."
    η::T3
    "current iteration."
    iteration::Int
    "whether the previous iteration was accepted."
    isaccept::Bool
end

AbstractMCMC.getparams(state::RobustAdaptiveMetropolis2State) = state.x
function AbstractMCMC.setparams!!(state::RobustAdaptiveMetropolis2State, x)
    return RobustAdaptiveMetropolis2State(
        x, state.logprob, state.S, state.logα, state.η, state.iteration, state.isaccept
    )
end

function ram_step_inner(
    rng::Random.AbstractRNG,
    model::AbstractMCMC.LogDensityModel,
    sampler::RobustAdaptiveMetropolis2,
    state::RobustAdaptiveMetropolis2State,
)
    # This is the initial state.
    f = model.logdensity
    d = LogDensityProblems.dimension(f)

    # Sample the proposal.
    x = state.x
    U = randn(rng, eltype(x), d)
    x_new = muladd(state.S, U, x)

    # Compute the acceptance probability.
    lp = state.logprob
    lp_new = LogDensityProblems.logdensity(f, x_new)
    # Technically, the `min` here is unnecessary for sampling according to `min(..., 1)`.
    # However, `ram_adapt` assumes that `logα` actually represents the log acceptance probability
    # and is thus bounded at 0. Moreover, users might be interested in inspecting the average
    # acceptance rate to check that the algorithm achieves something similar to the target rate.
    # Hence, it's a bit more convenient for the user if we just perform the `min` here
    # so they can just take an average of (`exp` of) the `logα` values.
    logα = min(lp_new - lp, zero(lp))
    isaccept = Random.randexp(rng) > -logα

    return x_new, lp_new, U, logα, isaccept
end

function ram_adapt(
    sampler::RobustAdaptiveMetropolis2,
    state::RobustAdaptiveMetropolis2State,
    logα::Real,
    U::AbstractVector,
)
    Δα = exp(logα) - sampler.α
    S = state.S
    # TODO: Make this configurable by defining a more general path.
    η = state.iteration^(-sampler.γ)
    ΔS = sqrt(η * abs(Δα)) * S * U / LinearAlgebra.norm(U)
    # TODO: Maybe do in-place and then have the user extract it with a callback if they really want it.
    S_new = if sign(Δα) == 1
        # One rank update.
        LinearAlgebra.lowrankupdate(LinearAlgebra.Cholesky(S.data, :L, 0), ΔS).L
    else
        # One rank downdate.
        LinearAlgebra.lowrankdowndate(LinearAlgebra.Cholesky(S.data, :L, 0), ΔS).L
    end
    return S_new, η
end

function AbstractMCMC.step(
    rng::Random.AbstractRNG,
    model::AbstractMCMC.LogDensityModel,
    sampler::RobustAdaptiveMetropolis2;
    initial_params=nothing,
    kwargs...,
)
    # This is the initial state.
    f = model.logdensity
    d = LogDensityProblems.dimension(f)

    # Initial parameter state.
    T = if initial_params === nothing
        eltype(sampler.γ)
    else
        Base.promote_type(eltype(sampler.γ), eltype(initial_params))
    end
    x = if initial_params === nothing
        randn(rng, T, d)
    else
        convert(AbstractVector{T}, initial_params)
    end
    # Initialize the Cholesky factor of the covariance matrix.
    S_data = if sampler.S === nothing
        LinearAlgebra.diagm(0 => ones(T, d))
    else
        # Check the dimensionality of the provided `S`.
        if size(sampler.S) != (d, d)
            throw(ArgumentError("The provided `S` has the wrong dimensionality."))
        end
        convert(AbstractMatrix{T}, sampler.S)
    end
    S = LinearAlgebra.LowerTriangular(S_data)

    # Construct the initial state.
    lp = LogDensityProblems.logdensity(f, x)
    state = RobustAdaptiveMetropolis2State(x, lp, S, zero(T), 0, 1, true)

    return AdvancedMH.Transition(x, lp, true), state
end

function AbstractMCMC.step(
    rng::Random.AbstractRNG,
    model::AbstractMCMC.LogDensityModel,
    sampler::RobustAdaptiveMetropolis2,
    state::RobustAdaptiveMetropolis2State;
    kwargs...,
)
    # Take the inner step.
    x_new, lp_new, U, logα, isaccept = ram_step_inner(rng, model, sampler, state)
    # Adapt the proposal.
    S_new, η = ram_adapt(sampler, state, logα, U)
    # Check that `S_new` has eigenvalues in the desired range.
    if !valid_eigenvalues(
        S_new, sampler.eigenvalue_lower_bound, sampler.eigenvalue_upper_bound
    )
        # In this case, we just keep the old `S` (p. 13 in Vihola, 2012).
        S_new = state.S
    end

    # Update state.
    state_new = RobustAdaptiveMetropolis2State(
        isaccept ? x_new : state.x,
        isaccept ? lp_new : state.logprob,
        S_new,
        logα,
        η,
        state.iteration + 1,
        isaccept,
    )
    return AdvancedMH.Transition(state_new.x, state_new.logprob, state_new.isaccept), state_new
end

function valid_eigenvalues(S, lower_bound, upper_bound)
    # Short-circuit if the bounds are the default.
    (lower_bound == 0 && upper_bound == Inf) && return true
    # Note that this is just the diagonal when `S` is triangular.
    eigenvals = LinearAlgebra.eigvals(S)
    return all(x -> lower_bound <= x <= upper_bound, eigenvals)
end

isgibbscomponent(::RobustAdaptiveMetropolis2) = true


AbstractMCMC.getparams(state::RobustAdaptiveMetropolis2State) = state.x
Turing.Inference.getparams(model, state::RobustAdaptiveMetropolis2State) = state.x
function AbstractMCMC.setparams!!(state::RobustAdaptiveMetropolis2State, θ)
    return RobustAdaptiveMetropolis2State(
        θ,
        state.logprob,
        state.S,
        state.logα,
        state.η,
        state.iteration,
        state.isaccept
    )
end