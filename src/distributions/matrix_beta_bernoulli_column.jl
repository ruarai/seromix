struct MatrixBetaBernoulliColumn{T} <: AbstractMatrixBernoulli
    alpha::T
    beta::T
    i::Int
    j::Int
end

function Distributions._rand!(rng::AbstractRNG, d::MatrixBetaBernoulliColumn, Y::AbstractMatrix)
    for j in 1:d.j
        p = rand(rng, Beta(d.alpha, d.beta)) # TODO account for births
        for i in 1:d.i
            Y[i,j] = rand(rng, Bernoulli(p))
        end
    end
end

function Distributions._logpdf(d::MatrixBetaBernoulliColumn, x::AbstractMatrix{<:Bool})
    logprob = 0.0

    for j in 1:d.j
        n_true = sum(view(x, :, j))
        n_total = d.i # TODO account for births

        logprob += logbeta(n_true + d.alpha, n_total - n_true + d.beta) - logbeta(d.alpha, d.beta)
    end

    return logprob
end



