
library(targets)
library(tarchetypes)
library(quarto)

suppressMessages(tar_source())

tar_option_set(packages = c("tidyverse", "arrow", "tidybayes", "bayesplot"))


sim_studies <- tibble(run_name = str_c("sim_study_", c("simple_1", "hanam_2018_1", "hanam_2018_2")))
data_studies <- tibble(run_name = c("hanam_2018"))


list(
  tar_target(hanam_data, read_hanam_data()),
  
  tar_target(hanam_data_file, save_hdf5(hanam_data, "runs/hanam_2018/model_data.hdf5"), format = "file"),
  tar_target(hanam_data_plots, plot_model_data(hanam_data, "hanam_2018")),
  
  tar_target(fit_model_jl, "src/fit_model.jl", format = "file"),
  
  tar_map(
    values = sim_studies,
    tar_target(run_dir, str_c("runs/", run_name, "/")),
    
    tar_target(model_data_file, str_c(run_dir, "model_data.hdf5"), format = "file"),
    tar_target(model_data, read_model_data(model_data_file)),
    tar_target(sim_study_plots, plot_sim_study(run_name, model_data)),
    
    tar_target(chain_file, call_fit_chain(run_name, run_dir, model_data_file, fit_model_jl), format = "file"),
    
    tar_target(chain_report, render_quarto(run_name, run_dir, model_data_file, chain_file,  "R/sim_study_chain_report.qmd"))
  ),
  
  tar_map(
    values = data_studies,
    tar_target(run_dir, str_c("runs/", run_name, "/")),
    
    tar_target(model_data_file, str_c(run_dir, "model_data.hdf5"), format = "file"),
    tar_target(chain_file, call_fit_chain(run_name, run_dir, model_data_file, fit_model_jl), format = "file"),
    
    tar_target(chain_report, render_quarto(run_name, run_dir, model_data_file, chain_file, "R/data_study_chain_report.qmd"))
  )

)