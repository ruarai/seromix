struct MatrixBernoulli{T} <: DiscreteMatrixDistribution
    p::T
    i::Int
    j::Int
end

Base.size(d::MatrixBernoulli) = (d.i, d.j)



function Distributions.rand(rng::AbstractRNG, d::MatrixBernoulli)
    Y = Matrix{Bool}(undef, size(d))

    for i in 1:d.i, j in 1:d.j
        Y[i,j] = rand(rng, Bernoulli(d.p))
    end

    return Y
end


function Distributions._rand!(rng::AbstractRNG, d::MatrixBernoulli, Y::AbstractMatrix)
    for i in 1:d.i, j in 1:d.j
        Y[i,j] = rand(rng, Bernoulli(d.p))
    end
end

function Distributions._logpdf(d::MatrixBernoulli, x::AbstractMatrix{<:Bool})
    pdf = 0.0

    @inbounds for j in 1:d.j
        @simd for i in 1:d.i
            pdf += d.p * x[i, j] + (1 - d.p[i, j]) * (1 - x[i, j])
        end
    end

    return log(pdf) 
end



