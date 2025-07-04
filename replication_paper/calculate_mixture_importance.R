

library(targets)
library(tidyverse)
library(loo)

tar_source()


chain <- tar_read(chain_model_comparison_hanam_2018_age_1) |> 
  filter(iteration > 25000)

model_data <- tar_read(hanam_2018_age)
inf_prior <- list(name = "BetaBernoulli", alpha = 1.0, beta = 1.0)

chain <- chain %>%
  filter(iteration > n_warmup)

fixed_params <- list(tau = 1e-10)

logp <- get_julia_function("pointwise_likelihood")(
  chain, model_data, "kucharski", inf_prior,
  fixed_params = fixed_params
)



elpd <- get_elpd(chain, model_data, inf_prior, "age_effect")

sum(elpd)

