
library(targets)
library(tarchetypes)
library(quarto)

tar_source("R/")

tar_option_set(packages = c("tidyverse", "arrow", "tidybayes", "bayesplot"))


sim_studies <- tribble(
  ~run_name,
  "simple_1",
  "hanam_2018_1"
)

list(
  tar_target(hanam_data, read_hanam_data()),
  
  tar_target(save_hanam_data, save_hdf5(hanam_data, "runs/hanam_2018/model_data.hdf5")),
  tar_target(hanam_data_plots, plot_model_data(hanam_data, "hanam_2018")),
  
  tar_target(fit_model_jl, "src/fit_model.jl", format = "file"),
  
  tar_map(
    values = sim_studies,
    tar_target(run_dir, str_c("runs/sim_study_", run_name, "/")),
    tar_target(model_data_file, str_c(run_dir, "model_data.hdf5"), format = "file"),
    tar_target(model_data, read_model_data(model_data_file)),
    tar_target(sim_study_plots, plot_sim_study(run_name, model_data)),
    
    tar_target(chain_file, call_fit_chain(run_name, run_dir, model_data_file, fit_model_jl), format = "file"),
    
    tar_target(quarto_report_filename, str_c("chain_report_", run_name, ".pdf")),
    tar_target(
      chain_report,
      {      
        file.remove(str_c(run_dir, "chain_report.pdf"))
        quarto::quarto_render(
          "reports/chain_report.qmd",
          output_file = str_c(run_name, ".pdf"),
          execute_params = list(model_data_file = model_data_file, chain_file = chain_file),
          quiet = TRUE
        )
        file.copy(str_c(run_name, ".pdf"), str_c(run_dir, "chain_report.pdf"))
        file.remove(str_c(run_name, ".pdf"))
      }

    )
  )
  

)