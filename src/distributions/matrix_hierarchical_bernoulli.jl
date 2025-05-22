
struct MatrixHierarchicalBernoulli{T} <: DiscreteMatrixDistribution
    row_means::AbstractVector{T}
    column_means::AbstractVector{T}

    i::Int
    j::Int
end

Base.size(d::MatrixHierarchicalBernoulli) = (d.i, d.j)

function Distributions.rand(rng::AbstractRNG, d::MatrixHierarchicalBernoulli)
    Y = Matrix{Bool}(undef, size(d))
    Distributions._rand!(rng, d, Y)
    return Y
end

function Distributions._rand!(rng::AbstractRNG, d::MatrixHierarchicalBernoulli, Y::AbstractMatrix)
    for i in 1:d.i, j in 1:d.j
        p = logistic(d.row_means[i] + d.column_means[j])

        Y[i,j] = rand(rng, Bernoulli(p))
    end
end

function Distributions._logpdf(d::MatrixHierarchicalBernoulli, x::AbstractMatrix{<:Bool})
    lpdf = 0.0

    for i in 1:d.i, j in 1:d.j
        p = logistic(d.row_means[i] + d.column_means[j])

        lpdf += logpdf(Bernoulli(p), x[i, j])
    end

    return lpdf
end



