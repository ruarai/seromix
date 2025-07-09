
library(targets)
library(tarchetypes)
library(crew)

suppressMessages(tar_source())

tar_option_set(
  controller = crew_controller_local(workers = 16)
)

default_n_iterations <- 75000
default_n_warmup <- 50000
default_n_chain <- 4


source("_targets_data_studies.R")
source("_targets_sim_studies.R")


c(targets_studies, targets_sim_studies)