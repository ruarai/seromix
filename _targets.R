
library(targets)
library(tarchetypes)
library(quarto)

suppressMessages(tar_source())

tar_option_set(packages = c("tidyverse", "arrow"))


list(
  # Load the original HaNam (2018) dataset
  tar_target(hanam_data, read_hanam_data()),
  tar_target(hanam_data_file, save_hdf5(hanam_data, "runs/hanam_2018/model_data.hdf5"), format = "file"),
  tar_target(hanam_data_plots, plot_model_data(hanam_data, "hanam_2018", plot_individuals = FALSE)),
  
  # Load the HaNam (2018) dataset with inferred age
  tar_target(hanam_data_age, read_hanam_data(use_inferred_age = TRUE)),
  tar_target(hanam_data_age_file, save_hdf5(hanam_data_age, "runs/hanam_2018_age/model_data.hdf5"), format = "file"),
  tar_target(hanam_data_age_plots, plot_model_data(hanam_data_age, "hanam_2018_age", plot_individuals = FALSE)),
  
  # Load the Fluscape (2009) dataset with neutralisation assay
  tar_target(fluscape_data_neuts, read_fluscape_data_neuts()),
  tar_target(fluscape_data_neuts_file, save_hdf5(fluscape_data_neuts, "runs/fluscape_2009_neuts/model_data.hdf5"), format = "file"),
  tar_target(fluscape_data_neuts_plots, plot_model_data(fluscape_data_neuts, "fluscape_2009_neuts")),
  
  # Load the Fluscape (2009) dataset with HI assay
  tar_target(fluscape_data_HI, read_fluscape_data_HI()),
  tar_target(fluscape_data_HI_file, save_hdf5(fluscape_data_HI, "runs/fluscape_2009_HI/model_data.hdf5"), format = "file"),
  tar_target(fluscape_data_HI_plots, plot_model_data(fluscape_data_HI, "fluscape_2009_HI"))
  
  
  
  
)