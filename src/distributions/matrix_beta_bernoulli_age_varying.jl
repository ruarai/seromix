struct MatrixBetaBernoulliAgeVarying{T} <: AbstractMatrixBernoulli
    alpha::T
    beta::T

    n_t_steps::Int
    n_subjects::Int

    subject_birth_ix::Vector{Int}
    n_alive_age::Vector{Int}
    max_age::Int
end



function MatrixBetaBernoulliAgeVarying(alpha, beta, sp::StaticModelParameters)
    max_age = maximum(sp.n_t_steps .- sp.subject_birth_ix)
    n_alive_age = zeros(Int, max_age)

    for age in 1:max_age
        n_alive = 0
        for ix_subject in 1:sp.n_subjects
            ix_t = sp.subject_birth_ix[ix_subject] + age
            if ix_t >= 1 && ix_t <= sp.n_t_steps
                n_alive += 1
            end
        end
        n_alive_age[age] = n_alive
    end
    
    return MatrixBetaBernoulliAgeVarying(alpha, beta, sp.n_t_steps, sp.n_subjects, sp.subject_birth_ix, n_alive_age, max_age)
end

function Distributions._rand!(rng::AbstractRNG, d::MatrixBetaBernoulliAgeVarying, Y::AbstractMatrix)
    Y .= false

    for age in 1:d.max_age
        p = rand(rng, Beta(d.alpha, d.beta))
        for ix_subject in 1:d.n_subjects
            ix_t = d.subject_birth_ix[ix_subject] + age
            if ix_t >= 1 && ix_t <= d.n_t_steps
                Y[ix_t, ix_subject] = rand(rng, Bernoulli(p))
            end
        end
    end
end

function Distributions._logpdf(d::MatrixBetaBernoulliAgeVarying, x::AbstractMatrix{<:Bool})
    logprob = 0.0

    for age in 1:d.max_age
        n_inf = 0
        for ix_subject in 1:d.n_subjects
            ix_t = d.subject_birth_ix[ix_subject] + age
            if ix_t >= 1 && ix_t <= d.n_t_steps
                n_inf += x[ix_t, ix_subject]
            end
        end
        n_alive = d.n_alive_age[age]

        logprob += logbeta(n_inf + d.alpha, n_alive - n_inf + d.beta) - logbeta(d.alpha, d.beta)
    end

    return logprob
end



