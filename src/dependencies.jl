using Pkg

if dirname(Base.active_project()) != pwd()
    Pkg.activate(".")
end


if VERSION != v"1.11.4"
    println("Julia version has changed (now $VERSION). Is this correct?")
end

using Distributions
using Turing
using AbstractMCMC
using DynamicPPL

# using Plots, StatsPlots

using LinearAlgebra
using StatsBase

using DataFrames
using QuackIO
# using Parquet

using Random
using LogDensityProblems

using AbstractPPL

using JLD2



include("antibody_kinetics.jl")

include("inference_model.jl")


import Turing.Inference: isgibbscomponent


include("mh_infection_sampler.jl")
include("mh_parameter_sampler.jl")

include("helper_functions.jl")
include("sampling_functions.jl")

include("distributions/titre_arraynormal.jl")
include("distributions/matrix_bernoulli.jl")


const const_titre_min = 0.0
const const_titre_max = 8.0

# function DynamicPPL.unflatten(vi::DynamicPPL.VarInfo, spl::AbstractMCMC.AbstractSampler, x::AbstractVector)
#     md = DynamicPPL.unflatten(vi.metadata, spl, x)

#     return DynamicPPL.VarInfo(
#         md,
#         Base.RefValue{DynamicPPL.float_type_with_fallback(eltype(x))}(DynamicPPL.getlogp(vi)),
#         Ref(DynamicPPL.get_num_produce(vi)),
#     )
# end