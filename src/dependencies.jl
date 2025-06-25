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
using StaticArrays
using ProgressMeter
using LoopVectorization

# using Plots, StatsPlots

using LinearAlgebra
using StatsBase
using StatsFuns
using SpecialFunctions

using DataFrames
using QuackIO

using Random
using LogDensityProblems

using AbstractPPL
using JLD2

using Printf

using SliceSampling


include("likelihood_model/make_waning_model.jl")
include("likelihood_model/antibody_kinetics_kucharski.jl")
include("likelihood_model/waning_model_kucharski.jl")


import Turing.Inference: isgibbscomponent
include("inference/mh_parameter_sampler.jl")
include("inference/mh_infection_sampler.jl")
include("inference/sampling_functions.jl")
include("inference/initial_params.jl")
include("proposals/infection_proposals.jl")

include("helper_functions.jl")

include("distributions/titre_arraynormal.jl")
include("distributions/matrix_bernoulli.jl")
include("distributions/matrix_beta_bernoulli.jl")
include("distributions/matrix_beta_bernoulli_column.jl")
include("distributions/matrix_beta_bernoulli_row.jl")


const const_titre_min = 0.0
const const_titre_max = 8.0