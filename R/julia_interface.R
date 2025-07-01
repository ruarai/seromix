
get_julia_fit_model <- function() {
  require(JuliaCall)
  
  Sys.setenv("JULIA_NUM_THREADS" = 8)
  julia <- julia_setup(JULIA_HOME = "/home/forecast/.asdf/shims/")
  
  julia_source("src/R_interface/fit_model.jl")
  
  julia_function("fit_model")
}
