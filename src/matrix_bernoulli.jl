struct MatrixBernoulli{T} <: DiscreteMatrixDistribution
    p::Matrix{T}
end

Base.size(d::MatrixBernoulli) = size(d.p)



function Distributions.rand(rng::AbstractRNG, d::MatrixBernoulli)
    Y = Matrix{Bool}(undef, size(d.p))

    for i in 1:size(d.p, 1), j in 1:size(d.p, 2)
        Y[i,j] = rand(rng, Bernoulli(d.p[i,j]))
    end

    return Y
end


function Distributions._rand!(rng::AbstractRNG, d::MatrixBernoulli, Y::AbstractMatrix)
    for i in 1:size(d.p, 1), j in 1:size(d.p, 2)
        Y[i,j] = rand(rng, Bernoulli(d.p[i,j]))
    end
end


function Distributions._rand!(rng::AbstractRNG, d::MatrixBernoulli, Y::AbstractMatrix)
    for i in 1:size(d.p, 1), j in 1:size(d.p, 2)
        Y[i,j] = rand(rng, Bernoulli(d.p[i,j]))
    end
end

function Distributions._logpdf(d::MatrixBernoulli, x::AbstractMatrix{<:Bool})
    pdf = 0.0

    @inbounds for j in 1:size(d.p, 2)
        @simd for i in 1:size(d.p, 1)
            pdf += d.p[i, j] * x[i, j] + (1 - d.p[i, j]) * (1 - x[i, j])
        end
    end

    return log(pdf) 
end



