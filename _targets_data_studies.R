

data_preprocessing <- list(
  # Load the original HaNam (2018) dataset
  tar_target(hanam_2018, read_hanam_data()),
  tar_target(hanam_2018_file, save_hdf5(hanam_2018, "runs/hanam_2018/model_data.hdf5"), format = "file"),
  tar_target(hanam_2018_plots, plot_model_data(hanam_2018, "hanam_2018", plot_individuals = FALSE)),
  
  # Load the HaNam (2018) dataset with inferred age
  tar_target(hanam_2018_age, read_hanam_data(use_inferred_age = TRUE)),
  tar_target(hanam_2018_age_file, save_hdf5(hanam_2018_age, "runs/hanam_2018_age/model_data.hdf5"), format = "file"),
  tar_target(hanam_2018_age_plots, plot_model_data(hanam_2018_age, "hanam_2018_age", plot_individuals = FALSE)),
  
  # Load the Fluscape (2009) dataset with neutralisation assay
  tar_target(fluscape_2009_neuts, read_fluscape_data_neuts()),
  tar_target(fluscape_2009_neuts_file, save_hdf5(fluscape_2009_neuts, "runs/fluscape_2009_neuts/model_data.hdf5"), format = "file"),
  tar_target(fluscape_2009_neuts_plots, plot_model_data(fluscape_2009_neuts, "fluscape_2009_neuts")),
  
  # Load the Fluscape (2009) dataset with HI assay
  tar_target(fluscape_2009_HI, read_fluscape_data_HI()),
  tar_target(fluscape_2009_HI_file, save_hdf5(fluscape_2009_HI, "runs/fluscape_2009_HI/model_data.hdf5"), format = "file"),
  tar_target(fluscape_2009_HI_plots, plot_model_data(fluscape_2009_HI, "fluscape_2009_HI"))
)


all_infection_priors <- list(
  matrix_bernoulli_10, matrix_bernoulli_30, matrix_bernoulli_50,
  matrix_beta_bernoulli_1_1, matrix_beta_bernoulli_1_1_tv, matrix_beta_bernoulli_1_1_sv,
  matrix_beta_bernoulli_2.5_8, matrix_beta_bernoulli_2.5_8_tv, matrix_beta_bernoulli_2.5_8_sv)

data_runs <- bind_rows(
  # Compare prior/proposal choices on Ha Nam (2018) study:
  expand_grid(
    exp_group = "prior_proposal",
    
    run_name = "hanam_2018",
    proposal_name = c("uncorrected", "corrected"),
    infection_prior = list(matrix_bernoulli_50, matrix_beta_bernoulli_1_1),
    
    sampler_name = "original",
    use_corrected_titre = FALSE,
    initial_params_name = "kucharski_data_study",
    
    # Filler variables to ensure column exists
    mixture_importance_sampling = FALSE,
    fixed_params = list(NULL)
  ),
  
  # Compare prior/proposal choices on fluscape studies:
  expand_grid(
    exp_group = "prior_proposal_fluscape",

    run_name = c("fluscape_2009_neuts", "fluscape_2009_HI"),
    fixed_params = list(list(omega = 0.5, mu_short = 1e-10, sigma_short = 1e-10)),
    proposal_name = c("uncorrected", "corrected"),
    infection_prior = list(matrix_bernoulli_50, matrix_beta_bernoulli_1_1),
    initial_params_name = "kucharski_data_study_fluscape",
    
    sampler_name = "original",
    use_corrected_titre = FALSE
  ),
  
  # Compare effect of age inclusion/exclusion
  expand_grid(
    exp_group = "age_inclusion",
    
    run_name = c("hanam_2018", "hanam_2018_age"),
    infection_prior = list(matrix_beta_bernoulli_1_1),
    initial_params_name = "kucharski_data_study"
  ),
  
  # Compare effect of titre correction:
  expand_grid(
    exp_group = "titre_correction",
    
    run_name = "hanam_2018",
    infection_prior = list(matrix_beta_bernoulli_1_1),
    initial_params_name = "kucharski_data_study",
    use_corrected_titre = c(TRUE, FALSE)
  ),
  
  # Compare priors:
  expand_grid(
    exp_group = "prior_comparison",
    
    run_name = "hanam_2018",
    infection_prior = all_infection_priors,
    initial_params_name = "kucharski_data_study"
  ),
  
  # Initial params
  expand_grid(
    exp_group = "initial_params",
    
    run_name = "hanam_2018",
    infection_prior = list(matrix_beta_bernoulli_1_1),
    initial_params_name = c("kucharski_data_study", "broad")
  ) |> 
    expand_grid(
      tribble(
        ~sampler_name, ~n_iterations, ~n_warmup,
        "default", default_n_iterations, default_n_warmup,
        "slice_sampler", 20000, 10000 
      )
    ),
  
  # Baseline MixIS score
  expand_grid(
    exp_group = "baseline_mixis",
    
    run_name = c("hanam_2018", "hanam_2018_age"),
    infection_prior = list(matrix_beta_bernoulli_1_1),
    initial_params_name = "kucharski_data_study",
    mixture_importance_sampling = TRUE,
    n_iterations = 200000, 
    n_warmup = 150000,
    n_chain = 4
  ),
  
  # Compare sampler
  tribble(
    ~sampler_name, ~n_iterations, ~n_warmup,
    "default", default_n_iterations, default_n_warmup,
    "slice_sampler", 20000, 10000 
  ) |> 
    mutate(
      exp_group = "sampler_comparison",
      run_name = "hanam_2018",
      infection_prior = list(matrix_beta_bernoulli_1_1),
      initial_params_name = "kucharski_data_study",
      n_chain = 4
    ),
  
  # Model comparison
  tribble(
    ~exp_name, ~initial_params_name, ~turing_model_name, ~fixed_params,
    "without_seniority", "age_effect", "age_effect", list(tau = 0.0, beta = 0.0, intercept = 0.0),
    "kucharski", "age_effect", "age_effect", list(beta = 0.0, intercept = 0.0),
    "with_age_effect", "age_effect", "age_effect", list(intercept = 0.0),
    "with_intercept", "age_effect", "age_effect", list(beta = 0.0),
    "with_both", "age_effect", "age_effect", NULL,
  ) |> 
    mutate(exp_group = "model_comparison",
           run_name = "hanam_2018_age",
           infection_prior = list(matrix_beta_bernoulli_1_1),
           n_iterations = 200000, 
           n_warmup = 150000,
           n_chain = 4
           ) |> 
    # Add mixture sampling to do model comparison
    expand_grid(mixture_importance_sampling = c(TRUE, FALSE))
  
) |>
  rowwise() |> 
  # Fill in defaults
  mutate(
    use_corrected_titre = replace_na(use_corrected_titre, TRUE),
    proposal_name = replace_na(proposal_name, "corrected"),
    turing_model_name = replace_na(turing_model_name, "kucharski"),
    mixture_importance_sampling = replace_na(mixture_importance_sampling, FALSE),
    sampler_name = replace_na(sampler_name, "default"),
    
    n_iterations = replace_na(n_iterations, default_n_iterations),
    n_warmup = replace_na(n_warmup, default_n_warmup),
    n_chain = replace_na(n_chain, default_n_chain),
  ) |>
  # Add a description of the prior
  mutate(prior_description = str_c(unlist(infection_prior), collapse = "_")) |>
  group_by(exp_group, run_name) |> 
  # Create a unique name and add in the run data (should probably be model_data)
  mutate(name = str_c(exp_group, "_", run_name, "_", row_number()),
         model_data = rlang::syms(run_name))


# Data which will be added to summarised outputs from each run
data_runs_meta <- data_runs |>
  select(name, exp_group, run_name, proposal_name, initial_params_name, 
         use_corrected_titre, prior_description, turing_model_name, 
         exp_name, mixture_importance_sampling, sampler_name)



data_chains <- tar_map(
  data_runs, names = name,
  tar_target(
    chain,
    get_julia_function("fit_model")(
      model_data,
      proposal_name = proposal_name,
      infection_prior = infection_prior,
      fixed_params = fixed_params,
      initial_params_name = initial_params_name,
      use_corrected_titre = use_corrected_titre,
      turing_model_name = turing_model_name,
      sampler_name = sampler_name,
      
      mixture_importance_sampling = mixture_importance_sampling,
      
      n_samples = as.integer(n_iterations),
      n_thinning = as.integer(round(n_iterations / 2000)),
      n_chain = as.integer(n_chain)
    ),
    garbage_collection = TRUE, memory = "transient", format = "parquet"
  ),
  
  tar_target(chain_subset, make_chain_subset(chain, model_data, name)),
  
  tar_target(chain_summary, summarise_chain(chain, n_warmup, model_data, add_name = name)),
  tar_target(chain_summary_singular, summarise_chain(chain, n_warmup, model_data, by_chain = FALSE, add_name = name)),
  
  tar_target(model_lp_mixis, get_lp_mixis(chain, n_warmup, model_data, turing_model_name, infection_prior, fixed_params, mixture_importance_sampling, name))
)


targets_studies <- list(
  data_preprocessing,
  data_chains,
  tar_combine(
    combined_summaries,
    data_chains[["chain_summary"]],
    command = dplyr::bind_rows(!!!.x) |> left_join(data_runs_meta, by = "name")
  ),
  tar_combine(
    combined_singular_summaries,
    data_chains[["chain_summary_singular"]],
    command = dplyr::bind_rows(!!!.x) |> left_join(data_runs_meta, by = "name")
  ),
  tar_combine(
    combined_chains,
    data_chains[["chain_subset"]],
    command = dplyr::bind_rows(!!!.x) |> left_join(data_runs_meta, by = "name")
  ),
  tar_combine(
    combined_lp_mixis,
    data_chains[["model_lp_mixis"]],
    command = dplyr::bind_rows(!!!.x) |> left_join(data_runs_meta, by = "name")
  )
)

