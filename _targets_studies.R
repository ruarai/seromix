

# Infection prior definitions (see src/R_interface/interface.jl)
matrix_bernoulli_10 <- list(name = "Bernoulli", p = 0.1)
matrix_bernoulli_30 <- list(name = "Bernoulli", p = 0.3)
matrix_bernoulli_50 <- list(name = "Bernoulli", p = 0.5)
matrix_beta_bernoulli_1_1 <- list(name = "BetaBernoulli", alpha = 1.0, beta = 1.0)
matrix_beta_bernoulli_1_1_tv <- list(name = "BetaBernoulliTimeVarying", alpha = 1.0, beta = 1.0)
matrix_beta_bernoulli_1_1_sv <- list(name = "BetaBernoulliSubjectVarying", alpha = 1.0, beta = 1.0)

matrix_beta_bernoulli_2.5_8 <- list(name = "BetaBernoulli", alpha = 2.5, beta = 8.0)
matrix_beta_bernoulli_2.5_8_tv <- list(name = "BetaBernoulliTimeVarying", alpha = 2.5, beta = 8.0)
matrix_beta_bernoulli_2.5_8_sv <- list(name = "BetaBernoulliSubjectVarying", alpha = 2.5, beta = 8.0)

data_runs <- bind_rows(
  # Compare prior/proposal choices on Ha Nam (2018) study:
  expand_grid(
    exp_group = "prior_proposal",
    
    run_name = "hanam_2018",
    fixed_params = list(NULL),sampler_name = "default", # filler
    proposal_name = c("uncorrected", "corrected"),
    infection_prior = list(matrix_bernoulli_50, matrix_beta_bernoulli_1_1),
    initial_params_name = "kucharski_data_study"
  ),
  
  # Compare prior/proposal choices on fluscape studies:
  # expand_grid(
  #   exp_group = "prior_proposal_fluscape",
  #   
  #   run_name = c("fluscape_2009_neuts", "fluscape_2009_HI"),
  #   fixed_params = list(list(omega = 0.5, mu_short = 1e-10, sigma_short = 1e-10)),
  #   proposal_name = c("uncorrected", "corrected"),
  #   infection_prior = list(matrix_bernoulli_50, matrix_beta_bernoulli_1_1),
  #   initial_params_name = "kucharski_data_study_fluscape"
  # ),
  
  # Compare effect of titre correction:
  expand_grid(
    exp_group = "titre_correction",
    
    run_name = "hanam_2018",
    infection_prior = list(matrix_beta_bernoulli_1_1),
    initial_params_name = "kucharski_data_study",
    use_corrected_titre = c(TRUE, FALSE)
  ),
  
  
  # Compare effect of titre correction on fluscape studies:
  # expand_grid(
  #   exp_group = "titre_correction",
  #   
  #   run_name = c("fluscape_2009_neuts", "fluscape_2009_HI"),
  #   fixed_params = list(list(omega = 0.5, mu_short = 1e-10, sigma_short = 1e-10)),
  #   infection_prior = list(matrix_beta_bernoulli_1_1),
  #   initial_params_name = "kucharski_data_study_fluscape",
  #   use_corrected_titre = c(TRUE, FALSE)
  # ),
  
  # Compare priors:
  expand_grid(
    exp_group = "prior_comparison",
    
    run_name = "hanam_2018",
    infection_prior = list(
      matrix_bernoulli_10, matrix_bernoulli_30, matrix_bernoulli_50,
      matrix_beta_bernoulli_1_1, matrix_beta_bernoulli_1_1_tv, matrix_beta_bernoulli_1_1_sv,
      matrix_beta_bernoulli_2.5_8, matrix_beta_bernoulli_2.5_8_tv, matrix_beta_bernoulli_2.5_8_sv),
    initial_params_name = "kucharski_data_study"
  ),
  
  # Compare effect of initial conditions:
  expand_grid(
    exp_group = "initial_conditions",
    
    run_name = "hanam_2018",
    infection_prior = list(matrix_beta_bernoulli_1_1),
    initial_params_name = c("broad", "kucharski_data_study")
  ),
  
  # Compare effect of sampler:
  # expand_grid(
  #   exp_group = "sampler_choice",
  #   
  #   run_name = "hanam_2018",
  #   infection_prior = list(matrix_beta_bernoulli_1_1, matrix_beta_bernoulli_2.5_8),
  #   initial_params_name = "broad",
  #   sampler_name = c("default", "slice_sampler")
  # )
) |>
  rowwise() |> 
  # Fill in defaults
  mutate(
    use_corrected_titre = replace_na(use_corrected_titre, TRUE),
    proposal_name = replace_na(proposal_name, "corrected"), # TODO NULL might work fine here?
    sampler_name = replace_na(sampler_name, "default")
  ) |>
  # Add a description of the prior
  mutate(prior_description = str_c(unlist(infection_prior), collapse = "_")) |>
  group_by(exp_group, run_name) |> 
  # Create a unique name and add in the run data (should probably be model_data)
  mutate(name = str_c(exp_group, "_", run_name, "_", row_number()),
         model_data = rlang::syms(run_name))


# Data which will be added to summarised outputs from each run
data_runs_meta <- data_runs |>
  select(name, exp_group, run_name, proposal_name, initial_params_name, use_corrected_titre, prior_description)


