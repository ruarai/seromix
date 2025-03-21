
library(targets)
library(tarchetypes)

tar_source("R/")

tar_option_set(packages = c("tidyverse", "arrow", "tidybayes", "bayesplot"))

list(
  tar_target(hanam_data, read_hanam_data()),
  tar_target(save_hanam_data, save_hdf5(hanam_data, "runs/hanam_2018/model_data.hdf5")),
  
  tar_target(hanam_data_plots, plot_model_data(hanam_data, "hanam_2018")),
  
  
  tar_target(sim_study_1_data, "runs/sim_study_simple_1/model_data.hdf5", format = "file"),
  tar_target(sim_study_1_plots, plot_sim_study("simple_1", sim_study_1_data))
  
  # tar_target(sim_study_hanam_2018_1_plots, prepare_sim_study("hanam_2018_1", hanam_data)),
  
  
  # tar_target(sim_study_hanam_2018_2_data, "runs/sim_study_hanam_2018_2/obs.parquet", format = tar_file()),
  # tar_target(sim_study_hanam_2018_2_plots, prepare_sim_study("hanam_2018_2", hanam_data)),
)