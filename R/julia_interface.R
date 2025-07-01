
setup_julia_env <- function() {
  require(JuliaCall)
  
  Sys.setenv("JULIA_NUM_THREADS" = 8)
  julia <- julia_setup(JULIA_HOME = "/home/forecast/.asdf/shims/")
  
  julia_source("src/R_interface/interface.jl")
}

get_julia_function <- function(function_name) {
  setup_julia_env()
  julia_function(function_name)
}


