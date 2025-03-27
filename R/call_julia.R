

call_fit_chain <- function(run_name, run_dir, data_file, fit_model_jl) {
  print(data_file)
  print(fit_model_jl)
  
  sim_study_name <- str_c("sim_study_", run_name)
  
  system(str_c("~/.asdf/shims/julia --threads=8 src/fit_model.jl ", sim_study_name))
  
  return(str_c(run_dir, "/chain.parquet"))
}
