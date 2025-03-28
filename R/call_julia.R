

call_fit_chain <- function(run_name, run_dir, data_file, fit_model_jl) {
  print(data_file)
  print(fit_model_jl)
  
  system(str_c("~/.asdf/shims/julia --threads=8 src/fit_model.jl ", run_name))
  
  return(str_c(run_dir, "/chain.parquet"))
}
