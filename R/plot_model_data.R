
plot_model_data <- function(model_data, data_name) {
  observations_df <- model_data$observations_df
  
  
  plot_dir <- str_c("results/", data_name, "/")
  dir.create(plot_dir, showWarnings = FALSE)
  
  
  
  sample_dates <- observations_df %>%
    distinct(subject_id, year_sampled)
  
  strains_observed <- observations_df %>%
    distinct(subject_id, strain_year)
  
  ggplot() +
    
    geom_rect(aes(xmin = -Inf, xmax = year_of_birth, ymin = subject_id - 0.5, ymax = subject_id + 0.5),
              fill = "grey90",
              model_data$subject_birth_data) +
    
    geom_tile(aes(x = strain_year, y = subject_id),
               size = 0.5, fill = "red", alpha = 0.3,
               strains_observed) +
    geom_point(aes(x = year_sampled, y = subject_id),
               size = 0.5,
               sample_dates) +
    
    geom_rug(aes(x = year), tibble(year = model_data$modelled_years)) +
    
    theme_bw()
  
  ggsave(
    str_c(plot_dir, "sample_dates.pdf"),
    bg = "white",
    width = 7, height = 5
  )
  
  
  ggplot() +
    geom_point(aes(x = X, y = Y, colour = strain_year), 
               model_data$raw_strain_coords) +
    
    geom_path(aes(x = x, y = y, colour = strain_year),
              model_data$fit_strain_coords) +
    
    geom_point(aes(x = x, y = y, colour = strain_year),
               model_data$fit_strain_coords %>% filter(strain_year %in% model_data$modelled_years)) +
    
    theme_bw()
  
  
  ggsave(
    str_c(plot_dir, "antigenic_map_fit.pdf"),
    bg = "white",
    width = 7, height = 5
  )
  
  
  n_pages <- floor(max(observations_df$subject_id) / 5)
  pdf(str_c(plot_dir, "data_by_year_sampled.pdf"), width = 10, height = 6)
  for(i in 1:n_pages) {
    p <- observations_df %>%
      filter(subject_id %in% (i * 5):((i + 1) * 5 - 1)) %>% 
      ggplot() +
      geom_point(aes(x = strain_year, y = observed_titre),
                 size = 0.7) +
      
      facet_grid(rows = vars(year_sampled), cols = vars(subject_id)) +
      
      coord_cartesian(ylim = c(0, 9)) +
      
      theme_bw()
    
    plot(p)
  }
  dev.off()
  
  
  pdf(str_c(plot_dir, "data_by_strain_year.pdf"), width = 10, height = 6)
  for(i in 1:max(observations_df$subject_id)) {
    p <- observations_df %>%
      filter(subject_id == i) %>% 
      ggplot() +
      geom_point(aes(x = year_sampled, y = observed_titre),
                 size = 0.7) +
      
      facet_wrap(~strain_year, ncol = 5) +
      
      coord_cartesian(ylim = c(0, 9)) +
      
      theme_bw()
    
    plot(p)
  }
  dev.off()
}

