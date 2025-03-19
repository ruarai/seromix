
library(targets)
library(tarchetypes)

tar_source("R/")

tar_option_set(packages = c("tidyverse", "arrow", "tidybayes", "bayesplot"))

list(
  tar_target(hanam_data, read_hanam_data()),
  tar_target(hanam_data_plots, plot_model_data(hanam_data, "hanam_2018")),
  
  
  
  
)