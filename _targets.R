
library(targets)
library(tarchetypes)
library(quarto)

suppressMessages(tar_source())

tar_option_set(packages = c("tidyverse", "arrow"))


list(
  tar_target(hanam_data, read_hanam_data()),
  tar_target(hanam_data_file, save_hdf5(hanam_data, "runs/hanam_2018/model_data.hdf5"), format = "file"),
  tar_target(hanam_data_plots, plot_model_data(hanam_data, "hanam_2018", plot_individuals = TRUE)),
  
  tar_target(fluscape_data_neuts, read_fluscape_data_neuts()),
  tar_target(fluscape_data_neuts_file, save_hdf5(fluscape_data_neuts, "runs/fluscape_2009_neuts/model_data.hdf5"), format = "file"),
  tar_target(fluscape_data_neuts_plots, plot_model_data(fluscape_data_neuts, "fluscape_2009_neuts")),
  
  tar_target(fluscape_data_HI, read_fluscape_data_HI()),
  tar_target(fluscape_data_HI_file, save_hdf5(fluscape_data_HI, "runs/fluscape_2009_HI/model_data.hdf5"), format = "file"),
  tar_target(fluscape_data_HI_plots, plot_model_data(fluscape_data_HI, "fluscape_2009_HI"))
  
)