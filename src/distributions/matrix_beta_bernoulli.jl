


struct MatrixBetaBernoulli{T} <: AbstractMatrixBernoulli
    alpha::T
    beta::T
    n_t_steps::Int
    n_subjects::Int

    subject_birth_ix::Vector{Int}

    n_total::Int
end

function MatrixBetaBernoulli(alpha, beta, p::StaticModelParameters)
    n_total = sum(sum(p.subject_birth_ix .<= ix_t) for ix_t in 1:p.n_t_steps)
    
    return MatrixBetaBernoulli(alpha, beta, p.n_t_steps, p.n_subjects, p.subject_birth_ix, n_total)
end

function Distributions._rand!(rng::AbstractRNG, d::MatrixBetaBernoulli, Y::AbstractMatrix)
    p = rand(rng, Beta(d.alpha, d.beta))

    for ix_t in 1:d.n_t_steps, ix_subject in 1:d.n_subjects
        if ix_t >= d.subject_birth_ix[ix_subject]
            Y[ix_t, ix_subject] = rand(rng, Bernoulli(p))
        else
            Y[ix_t, ix_subject] = false
        end
    end
end

function Distributions._logpdf(d::MatrixBetaBernoulli, x::AbstractMatrix{<:Bool})
    n_true = sum(x)
    n_total = d.n_total

    return logbeta(n_true + d.alpha, n_total - n_true + d.beta) - logbeta(d.alpha, d.beta)
end



