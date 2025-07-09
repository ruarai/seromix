


include("dependencies.jl")


data_code = "hanam_2018_age"

run_dir = "runs/$(data_code)/"

model_data = load("$run_dir/model_data.hdf5")

p = read_model_parameters(model_data)

chain_name = "nonlinear_test"

chain_df = QuackIO.read_parquet(DataFrame, "$run_dir/chain_$chain_name.parquet")

chain_df_filt = filter(:iteration => iteration -> iteration > 15000, chain_df)

ppd_df = model_ppd(chain_df_filt, p, model; n_samples = 30)

save_draws(ppd_df, "$run_dir/ppd_$(chain_name).parquet")
