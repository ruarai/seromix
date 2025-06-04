

include("../dependencies.jl")
include("../inference/summarise_chain.jl")



model_data = load("runs/sim_study_hanam_2018_4/1/model_data.hdf5")

chain_files = list_files("runs/sim_study_hanam_2018_4", r"^chain_.*\.parquet$")

for i in eachindex(chain_files)
    summary = summarise_chain(chain_files[i], 50_000, model_data)


    write_table(replace(chain_files[i], ".parquet" => "_summary.csv"), summary.chain_summaries)
    write_table(replace(chain_files[i], ".parquet" => "_diagnostics.csv"), summary.diagnostics)
end


