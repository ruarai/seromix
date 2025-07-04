


sim_ar_runs <- expand_grid(
  exp_group = "attack_rate_prior",
  endemic_mean_ar = c(0.1, 0.3, 0.5),
  infection_prior = list(matrix_beta_bernoulli_1_1, matrix_beta_bernoulli_1_1_tv, matrix_beta_bernoulli_1_1_sv),
  drop_age = c(TRUE, FALSE)
  # infection_prior = list(
  #   matrix_beta_bernoulli_1_1, matrix_beta_bernoulli_1_1_tv, matrix_beta_bernoulli_1_1_sv,
  #   matrix_beta_bernoulli_2.5_8, matrix_beta_bernoulli_2.5_8_tv, matrix_beta_bernoulli_2.5_8_sv)
) |>
  rowwise() |> 
  mutate(prior_description = str_c(unlist(infection_prior), collapse = "_")) |> 
  ungroup() |> 
  mutate(name = str_c(exp_group, "_", endemic_mean_ar, "_", prior_description, "_", row_number()))

sim_ar_meta <- sim_ar_runs |> 
  select(name, exp_group, endemic_mean_ar, prior_description, drop_age)


sim_ar_chains <- tar_map(
  sim_ar_runs, names = name,
  
  tar_target(sim_model_data, get_julia_function("simulate_hanam_2018")(endemic_mean_ar, drop_age = drop_age)),
  
  tar_target(
    chain,
    get_julia_function("fit_model")(
      sim_model_data,
      infection_prior = infection_prior,
      initial_params_name = "kucharski_sim_study",
      n_samples = as.integer(default_n_iterations),
      n_thinning = as.integer(round(default_n_iterations / 2000)),
      n_chain = as.integer(default_n_chain)
    ),
    garbage_collection = TRUE, format = "parquet"
  ),
  tar_target(chain_subset, make_chain_subset(chain, sim_model_data, name)),
  
  tar_target(chain_summary, summarise_chain(chain, default_n_warmup, sim_model_data, add_name = name)),
  tar_target(chain_summary_singular, summarise_chain(chain, default_n_warmup, sim_model_data, by_chain = FALSE, add_name = name))
)

targets_sim_studies <- list(
  sim_ar_chains,
  tar_combine(
    combined_sim_ar_summaries,
    sim_ar_chains[["chain_summary"]],
    command = dplyr::bind_rows(!!!.x) |> left_join(sim_ar_meta, by = "name")
  ),
  tar_combine(
    combined_sim_ar_singular_summaries,
    sim_ar_chains[["chain_summary_singular"]],
    command = dplyr::bind_rows(!!!.x) |> left_join(sim_ar_meta, by = "name")
  )
)






