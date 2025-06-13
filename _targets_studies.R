

# Infection prior definitions (see src/R_interface/interface.jl)
matrix_bernoulli_50 <- list(name = "Bernoulli", p = 0.5)
matrix_beta_bernoulli_1_1 <- list(name = "BetaBernoulli", alpha = 1.0, beta = 1.0)
matrix_beta_bernoulli_1_1_tv <- list(name = "BetaBernoulliTimeVarying", alpha = 1.0, beta = 1.0)
matrix_beta_bernoulli_1_1_sv <- list(name = "BetaBernoulliSubjectVarying", alpha = 1.0, beta = 1.0)

base_hanam_2018 <- tibble(run_name = "hanam_2018", fixed_params = list(NULL))

data_runs <- bind_rows(
  # Compare prior/proposal choices on Ha Nam (2018) study:
  expand_grid(
    exp_group = "prior_proposal",
    
    run_name = "hanam_2018",
    fixed_params = list(NULL),
    proposal_name = c("uncorrected", "corrected"),
    infection_prior = list(matrix_bernoulli_50, matrix_beta_bernoulli_1_1,
                           matrix_beta_bernoulli_1_1_tv, matrix_beta_bernoulli_1_1_sv),
    initial_params_name = "kucharski_data_study",
    use_corrected_titre = TRUE
  ),
  
  # Compare prior/proposal choices on fluscape studies:
  expand_grid(
    exp_group = "prior_proposal_fluscape",
    
    run_name = c("fluscape_2009_neuts", "fluscape_2009_HI"),
    fixed_params = list(NULL),
    proposal_name = c("uncorrected", "corrected"),
    infection_prior = list(matrix_bernoulli_50, matrix_beta_bernoulli_1_1,
                           matrix_beta_bernoulli_1_1_tv, matrix_beta_bernoulli_1_1_sv),
    initial_params_name = "kucharski_data_study_fluscape",
    use_corrected_titre = TRUE
  ),
  
  # Compare effect of titre correction:
  expand_grid(
    exp_group = "titre_correction",
    
    run_name = "hanam_2018",
    fixed_params = list(NULL),
    proposal_name = "corrected",
    infection_prior = list(matrix_beta_bernoulli_1_1),
    initial_params_name = "kucharski_data_study",
    use_corrected_titre = c(TRUE, FALSE)
  )
  
  # Compare effect of initial conditions:
  
  
  # Try using slice sampling:
  
  
  
) %>%
  rowwise() %>%
  mutate(prior_description = str_c(unlist(infection_prior), collapse = "_")) %>%
  ungroup() %>%
  mutate(name = str_c(run_name, "_", row_number()),
         run_data = rlang::syms(run_name))



data_runs_meta <- data_runs %>%
  select(name, exp_group, run_name, proposal_name, initial_params_name, use_corrected_titre, prior_description)


