using Pkg

if dirname(Base.active_project()) != pwd()
    Pkg.activate(".")
end


if VERSION != v"1.11.3"
    println("Julia version has changed. Is this correct?")
end

using Distributions
using Turing

using Plots, StatsPlots

using LinearAlgebra
using StatsBase

using DataFrames
using Parquet

using Random
using LogDensityProblems



include("antibody_kinetics.jl")
include("inference_model.jl")
include("data_interface.jl")
# include("function_proposal.jl")