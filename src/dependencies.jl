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



# include("antibody_kinetics.jl")
include("antibody_kinetics_opt.jl")
include("inference_model.jl")
include("data_interface.jl")
# include("function_proposal.jl")
include("mh_sampler.jl")

# Expand grid via
# https://stackoverflow.com/a/67733908
function expand_grid(; iters...)
    var_names = collect(keys(iters))

    var_itr = [1:length(x) for x in values(iters)]
    var_ix = vcat([collect(x)' for x in Iterators.product(var_itr...)]...)
    out = DataFrame()
    for i = eachindex(var_names)
        out[:,var_names[i]] = collect(iters[i])[var_ix[:,i]]
    end
    return out
end