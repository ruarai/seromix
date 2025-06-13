struct MatrixBetaBernoulliRow{T} <: AbstractMatrixBernoulli
    alpha::T
    beta::T
    i::Int
    j::Int
end

function Distributions._rand!(rng::AbstractRNG, d::MatrixBetaBernoulliRow, Y::AbstractMatrix)
    for i in 1:d.i
        p = rand(rng, Beta(d.alpha, d.beta)) # TODO account for births
        for j in 1:d.j
            Y[i,j] = rand(rng, Bernoulli(p))
        end
    end
end

function Distributions._logpdf(d::MatrixBetaBernoulliRow, x::AbstractMatrix{<:Bool})
    logprob = 0.0

    for i in 1:d.i
        n_true = sum(view(x, i, :))
        n_total = d.j # TODO account for births

        logprob += logbeta(n_true + d.alpha, n_total - n_true + d.beta) - logbeta(d.alpha, d.beta)
    end

    return logprob
end



