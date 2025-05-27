


include("dependencies.jl")


data_code = "hanam_2018"

run_dir = "runs/$(data_code)/"

model_data = load("$run_dir/model_data.hdf5")

p = read_model_parameters(model_data)

chain_name = "infer_distances_1"

chain_df = QuackIO.read_parquet(DataFrame, "$run_dir/chain_$chain_name.parquet")

# chain.tau = 0.0

chain_filt = filter(:iteration => iteration -> iteration > 30_000, chain_df)



n_draws = 30
n_subjects = p.n_subjects

# ppd_df = ppd_kucharski(chain_filt, p, n_subjects, n_draws)
ppd_df = ppd_infer_distances(chain_filt, p, n_subjects, n_draws)

save_draws(ppd_df, "$run_dir/ppd_$(chain_name).parquet")

