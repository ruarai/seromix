struct MatrixBetaBernoulliSubjectVarying{T} <: AbstractMatrixBernoulli
    alpha::T
    beta::T
    n_t_steps::Int
    n_subjects::Int

    subject_birth_ix::Vector{Int}
end


function MatrixBetaBernoulliSubjectVarying(alpha, beta, p::FixedModelParameters)
    return MatrixBetaBernoulliSubjectVarying(alpha, beta, p.n_t_steps, p.n_subjects, p.subject_birth_ix)
end

function Distributions._rand!(rng::AbstractRNG, d::MatrixBetaBernoulliSubjectVarying, Y::AbstractMatrix)
    for ix_subject in 1:d.n_subjects
        p = rand(rng, Beta(d.alpha, d.beta))
        for ix_t in 1:d.n_t_steps
            if ix_t >= d.subject_birth_ix[ix_subject]
                Y[ix_t, ix_subject] = rand(rng, Bernoulli(p))
            else
                Y[ix_t, ix_subject] = false
            end
        end
    end
end

function Distributions._logpdf(d::MatrixBetaBernoulliSubjectVarying, x::AbstractMatrix{<:Bool})
    logprob = 0.0

    for ix_subject in 1:d.n_subjects
        ix_t_earliest = @inbounds max(1, d.subject_birth_ix[ix_subject])

        n_true = sum(view(x, ix_t_earliest:d.n_t_steps, ix_subject))
        n_total = d.n_t_steps - ix_t_earliest + 1

        logprob += logbeta(n_true + d.alpha, n_total - n_true + d.beta) - logbeta(d.alpha, d.beta)
    end

    return logprob
end



