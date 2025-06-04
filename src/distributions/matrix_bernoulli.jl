struct MatrixBernoulli{T} <: DiscreteMatrixDistribution
    p::T
    i::Int
    j::Int
end

Base.size(d::MatrixBernoulli) = (d.i, d.j)

Base.eltype(d::MatrixBernoulli) = Bool

function Distributions.rand(rng::AbstractRNG, d::MatrixBernoulli)
    Y = Matrix{Bool}(undef, size(d))
    Distributions._rand!(rng, d, Y)
    return Y
end

function Distributions._rand!(rng::AbstractRNG, d::MatrixBernoulli, Y::AbstractMatrix)
    for i in 1:d.i, j in 1:d.j
        Y[i,j] = rand(rng, Bernoulli(d.p))
    end
end

function Distributions._logpdf(d::MatrixBernoulli, x::AbstractMatrix{<:Bool})
    log_p = log(d.p)
    log_1_minus_p = log(1.0 - d.p)

    n_true = sum(x)
    n_false = (d.i * d.j) - n_true

    return n_true * log_p + n_false * log_1_minus_p
end
