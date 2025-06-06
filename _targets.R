
library(targets)
library(tarchetypes)
library(crew)

suppressMessages(tar_source())

tar_option_set(
  controller = crew_controller_local(workers = 8)
)

source("_targets_studies.R")

preprocessing <- list(
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

runs <- tar_map(
  data_runs, names = name,
  tar_target(
    chain,
    get_julia_fit_model()(
      run_data,
      proposal_name = proposal_name,
      infection_prior = infection_prior,
      fixed_params = fixed_params,
      initial_params_name = "kucharski_data_study",
      
      n_samples = as.integer(100000),
      n_thinning = as.integer(50),
      n_chain = as.integer(4)
    ),
    format = "parquet"
  ),
  tar_target(
    chain_summary,
    summarise_chain(chain, 70000, run_data) %>%
      mutate(name = name, 
             run_name = run_name, 
             proposal_name = proposal_name, 
             prior_description = prior_description)
  )
)

combine <- tar_combine(
  combined_summaries,
  runs[["chain_summary"]],
  command = dplyr::bind_rows(!!!.x)
)

list(
  preprocessing,
  runs,
  combine
)
