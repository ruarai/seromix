struct MatrixBetaBernoulliTimeVarying{T} <: AbstractMatrixBernoulli
    alpha::T
    beta::T

    n_t_steps::Int
    n_subjects::Int

    subject_birth_ix::Vector{Int}
    n_alive_t::Vector{Int}
end



function MatrixBetaBernoulliTimeVarying(alpha, beta, p::StaticModelParameters)
    n_alive_t = [sum(p.subject_birth_ix .<= ix_t) for ix_t in 1:p.n_t_steps]
    
    return MatrixBetaBernoulliTimeVarying(alpha, beta, p.n_t_steps, p.n_subjects, p.subject_birth_ix, n_alive_t)
end

function Distributions._rand!(rng::AbstractRNG, d::MatrixBetaBernoulliTimeVarying, Y::AbstractMatrix)
    for ix_t in 1:d.n_t_steps
        p = rand(rng, Beta(d.alpha, d.beta))

        for ix_subject in 1:d.n_subjects
            if ix_t >= d.subject_birth_ix[ix_subject]
                Y[ix_t, ix_subject] = rand(rng, Bernoulli(p))
            else
                Y[ix_t, ix_subject] = false
            end
        end
    end
end

function Distributions._logpdf(d::MatrixBetaBernoulliTimeVarying, x::AbstractMatrix{<:Bool})
    logprob = 0.0

    for ix_t in 1:d.n_t_steps
        n_true = sum(view(x, ix_t, :))
        n_total = @inbounds d.n_alive_t[ix_t]

        logprob += logbeta(n_true + d.alpha, n_total - n_true + d.beta) - logbeta(d.alpha, d.beta)
    end

    return logprob
end



