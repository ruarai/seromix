

library(targets)
library(tidyverse)
library(loo)

tar_source()


chain <- arrow::read_parquet("runs/hanam_2018/chain_loo_test.parquet") %>%
  filter(iteration > 25000)

model_data <- tar_read(hanam_2018)

logp <- get_julia_function("pointwise_likelihood_kucharski")(chain, model_data)

dim(logp)

library(matrixStats)


elpd_chains <- map(
  1:4,
  function(ix_chain) {
    chain_samples <- chain$chain == ix_chain
    l_common_mix <- rowLogSumExps(-logp[chain_samples, ])
    log_weights <- -logp[chain_samples, ] - l_common_mix
    elpd_mixis <- logSumExp(-l_common_mix) - rowLogSumExps(t(log_weights))
  }
)

elpd_mat <- do.call(cbind, elpd_chains)
chain_weights <- stacking_weights(elpd_mat)

chain$stacking_weights <- chain_weights[chain$chain]



ggdist::weighted_quantile(
  chain$sigma_long,
  probs = c(0.05, 0.5, 0.95),
  weights = chain$stacking_weights
)


