
abstract type AbstractMatrixBernoulli <: DiscreteMatrixDistribution end
Base.size(d::AbstractMatrixBernoulli) = (d.n_t_steps, d.n_subjects)
Base.eltype(d::AbstractMatrixBernoulli) = Bool


function Distributions.rand(rng::AbstractRNG, d::AbstractMatrixBernoulli)
    Y = Matrix{Bool}(undef, size(d))
    Distributions._rand!(rng, d, Y)
    return Y
end

# Necessary to get logjoint to work
function Distributions._logpdf(d::AbstractMatrixBernoulli, x::AbstractMatrix{<:Float64})
    return Distributions._logpdf(d, convert.(Bool, x))
end


struct MatrixBernoulli{T} <: AbstractMatrixBernoulli
    p::T
    n_t_steps::Int
    n_subjects::Int

    subject_birth_ix::Vector{Int}

    n_total::Int
end

function MatrixBernoulli(prob, p::FixedModelParameters)
    n_total = sum(sum(p.subject_birth_ix .<= ix_t) for ix_t in 1:p.n_t_steps)
    
    return MatrixBernoulli(prob, p.n_t_steps, p.n_subjects, p.subject_birth_ix, n_total)
end


function Distributions._rand!(rng::AbstractRNG, d::MatrixBernoulli, Y::AbstractMatrix)
    for ix_t in 1:d.n_t_steps, ix_subject in 1:d.n_subjects
        if ix_t >= d.subject_birth_ix[ix_subject]
            Y[ix_t, ix_subject] = rand(rng, Bernoulli(d.p))
        else
            Y[ix_t, ix_subject] = false
        end
    end
end

function Distributions._logpdf(d::MatrixBernoulli, x::AbstractMatrix{<:Bool})
    log_p = log(d.p)
    log_1_minus_p = log(1.0 - d.p)

    n_true = sum(x)
    n_false = d.n_total - n_true

    return n_true * log_p + n_false * log_1_minus_p
end
