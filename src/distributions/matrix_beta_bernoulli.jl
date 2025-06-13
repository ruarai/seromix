


struct MatrixBetaBernoulli{T} <: AbstractMatrixBernoulli
    alpha::T
    beta::T
    i::Int
    j::Int
end

function Distributions._rand!(rng::AbstractRNG, d::MatrixBetaBernoulli, Y::AbstractMatrix)
    p = rand(rng, Beta(d.alpha, d.beta))

    for i in 1:d.i, j in 1:d.j
        Y[i,j] = rand(rng, Bernoulli(p))
    end
end

function Distributions._logpdf(d::MatrixBetaBernoulli, x::AbstractMatrix{<:Bool})
    n_true = sum(x)
    n_total = d.i * d.j # TODO account for births

    return logbeta(n_true + d.alpha, n_total - n_true + d.beta) - logbeta(d.alpha, d.beta)
end



