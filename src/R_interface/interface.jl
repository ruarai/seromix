
include("../dependencies.jl")

include("fit_model.jl")
include("create_sim_study.jl")
include("sim_study_functions.jl")
include("pointwise_likelihood.jl")


function convert_model_data(model_data_R)
    return OrderedDict{String, Any}(String(k) => v for (k, v) in model_data_R)
end