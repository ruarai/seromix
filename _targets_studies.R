

# Infection prior definitions (see src/R_interface/interface.jl)
matrix_bernoulli_50 <- list(name = "Bernoulli", p = 0.5)
matrix_beta_bernoulli_1_1 <- list(name = "BetaBernoulli", alpha = 1.0, beta = 1.0)


data_runs <- bind_rows(
  tibble(run_name = "hanam_2018"),
  expand_grid(
    run_name = c("fluscape_2009_HI", "fluscape_2009_neuts"),
    fixed_params = list(list(omega = 0.5, mu_short = 1e-10, sigma_short = 1e-10))
  )
) %>% expand_grid(
  proposal_name = c("uncorrected", "corrected"),
  infection_prior = list(matrix_bernoulli_50, matrix_beta_bernoulli_1_1)
) %>%
  rowwise() %>%
  mutate(prior_description = str_c(unlist(infection_prior), collapse = "_")) %>%
  ungroup() %>%
  mutate(name = str_c(run_name, proposal_name, prior_description, sep = "_"),
         run_data = rlang::syms(run_name))




