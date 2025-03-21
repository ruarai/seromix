

plot_sim_study <- function(
  study_name,
  model_data
) {
  results_name <- str_c("sim_study_", study_name)
  run_dir <- str_c("runs/", results_name, "/")
  
  plot_model_data(model_data, results_name, plot_individuals = FALSE)
  
  
  infections <- model_data$infections
  
  
  sample_dates <- model_data$observations %>%
    distinct(ix_subject, year_observed)
  
  ggplot() +
    
    geom_rect(aes(xmin = -Inf, xmax = year_of_birth + 0.5, ymin = ix_subject - 0.5, ymax = ix_subject + 0.5),
              fill = "grey90",
              model_data$subject_birth_data) +
    
    geom_tile(aes(x = year, y = ix_subject),
              fill = "lightblue",
              infections) +
    
    geom_point(aes(x = year_observed, y = ix_subject),
               size = 0.5,
               sample_dates) +
    
    geom_rug(aes(x = year), tibble(year = model_data$modelled_years)) +
    
    theme_bw()
  
  
  ggsave(
    str_c(run_dir, "sample_dates_infections.pdf"),
    bg = "white",
    width = 7, height = 5
  )
}