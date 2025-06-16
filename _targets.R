
library(targets)
library(tarchetypes)
library(crew)

suppressMessages(tar_source())

tar_option_set(
  controller = crew_controller_local(workers = 3)
)

n_iterations <- 200000
n_warmup <- 150000
n_chain <- 8


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
      initial_params_name = initial_params_name,
      use_corrected_titre = use_corrected_titre,
      
      n_samples = as.integer(n_iterations),
      n_thinning = as.integer(round(n_iterations / 2000)),
      n_chain = as.integer(n_chain)
    ),
    garbage_collection = TRUE,
    format = "parquet"
  ),
  tar_target(
    chain_subset,
    chain |>
      select(-starts_with("infections")) |> 
      clean_chain() |>
      add_total_infections(run_data) |>
      mutate(name = name)
  ),
  tar_target(
    chain_summary,
    summarise_chain(chain, n_warmup, run_data, by_chain = TRUE) |>
      mutate(name = name)
  ),
  tar_target(
    chain_summary_singular,
    summarise_chain(chain, n_warmup, run_data, by_chain = FALSE) |>
      mutate(name = name)
  )
)


list(
  preprocessing,
  runs,
  tar_combine(
    combined_summaries,
    runs[["chain_summary"]],
    command = dplyr::bind_rows(!!!.x) |> left_join(data_runs_meta, by = "name")
  ),
  tar_combine(
    combined_singular_summaries,
    runs[["chain_summary_singular"]],
    command = dplyr::bind_rows(!!!.x) |> left_join(data_runs_meta, by = "name")
  ),
  tar_combine(
    combined_chains,
    runs[["chain_subset"]],
    command = dplyr::bind_rows(!!!.x) |> left_join(data_runs_meta, by = "name")
  )
)
