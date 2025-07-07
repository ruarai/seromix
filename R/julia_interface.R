

# Infection prior definitions (see src/R_interface/fit_model.jl)
matrix_bernoulli_10 <- list(name = "Bernoulli", p = 0.1)
matrix_bernoulli_30 <- list(name = "Bernoulli", p = 0.3)
matrix_bernoulli_50 <- list(name = "Bernoulli", p = 0.5)
matrix_beta_bernoulli_1_1 <- list(name = "BetaBernoulli", alpha = 1.0, beta = 1.0)
matrix_beta_bernoulli_1_1_tv <- list(name = "BetaBernoulliTimeVarying", alpha = 1.0, beta = 1.0)
matrix_beta_bernoulli_1_1_sv <- list(name = "BetaBernoulliSubjectVarying", alpha = 1.0, beta = 1.0)

matrix_beta_bernoulli_2.5_8 <- list(name = "BetaBernoulli", alpha = 2.5, beta = 8.0)
matrix_beta_bernoulli_2.5_8_tv <- list(name = "BetaBernoulliTimeVarying", alpha = 2.5, beta = 8.0)
matrix_beta_bernoulli_2.5_8_sv <- list(name = "BetaBernoulliSubjectVarying", alpha = 2.5, beta = 8.0)


setup_julia_env <- function() {
  require(JuliaCall)
  
  Sys.setenv("JULIA_NUM_THREADS" = 8)
  julia <- julia_setup(JULIA_HOME = "/home/forecast/.asdf/shims/")
  
  julia_source("src/R_interface/interface.jl")
}

# Reloads the Julia scripts and returns the given function
get_julia_function <- function(function_name) {
  setup_julia_env()
  julia_function(function_name)
}

