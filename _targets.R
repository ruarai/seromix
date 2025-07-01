
library(targets)
library(tarchetypes)
library(crew)

suppressMessages(tar_source())

tar_option_set(
  controller = crew_controller_local(workers = 32)
)

n_iterations <- 75000
n_warmup <- 50000
n_chain <- 4


source("_targets_studies.R")
source("_targets_sim_studies.R")



list(
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

