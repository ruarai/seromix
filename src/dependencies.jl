using Pkg

if dirname(Base.active_project()) != pwd()
    Pkg.activate(".")
end


if VERSION != v"1.11.4"
    println("Julia version has changed (now $VERSION). Is this correct?")
end

using Distributions
using Turing

# using Plots, StatsPlots

using LinearAlgebra
using StatsBase

using DataFrames
# using Parquet

using Random
using LogDensityProblems

using AbstractPPL

using JLD2



include("antibody_kinetics.jl")

include("inference_model.jl")
include("mh_sampler.jl")

include("helper_functions.jl")

# include("make_ppd.jl")
include("distributions/titre_arraynormal.jl")

# include("distributions/matrix_bernoulli.jl")


const const_titre_min = 0.0
const const_titre_max = 8.0