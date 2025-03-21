
plot_model_data <- function(model_data, run_name, plot_individuals = FALSE) {
  observations <- model_data$observations
  
  
  plot_dir <- str_c("runs/", run_name, "/")
  dir.create(plot_dir, showWarnings = FALSE, recursive = TRUE)
  
  
  
  sample_dates <- observations %>%
    distinct(ix_subject, year_observed)
  
  strains_observed <- observations %>%
    distinct(ix_subject, strain_year)
  
  ggplot() +
    
    geom_rect(aes(xmin = -Inf, xmax = year_of_birth, ymin = ix_subject - 0.5, ymax = ix_subject + 0.5),
              fill = "grey90",
              model_data$subject_birth_data) +
    
    geom_tile(aes(x = strain_year, y = ix_subject),
               fill = "red", alpha = 0.3,
               strains_observed) +
    geom_point(aes(x = year_observed, y = ix_subject),
               sample_dates) +
    
    geom_rug(aes(x = year), tibble(year = model_data$modelled_years)) +
    
    theme_bw()
  
  ggsave(
    str_c(plot_dir, "sample_dates.pdf"),
    bg = "white",
    width = 7, height = 5
  )
  
  if(!is.null(model_data$raw_strain_coords)) {  
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
  }
  

  
  if(plot_individuals) {
    n_pages <- floor(max(observations$ix_subject) / 5)
    pdf(str_c(plot_dir, "data_by_year_observed.pdf"), width = 10, height = 6)
    for(i in 1:n_pages) {
      p <- observations %>%
        filter(ix_subject %in% (i * 5):((i + 1) * 5 - 1)) %>% 
        ggplot() +
        geom_point(aes(x = strain_year, y = observed_titre),
                   size = 0.7) +
        
        facet_grid(rows = vars(year_observed), cols = vars(ix_subject)) +
        
        coord_cartesian(ylim = c(0, 9)) +
        
        theme_bw()
      
      plot(p)
    }
    dev.off()
    
    
    pdf(str_c(plot_dir, "data_by_strain_year.pdf"), width = 10, height = 6)
    for(i in 1:max(observations$ix_subject)) {
      p <- observations %>%
        filter(ix_subject == i) %>% 
        ggplot() +
        geom_point(aes(x = year_observed, y = observed_titre),
                   size = 0.7) +
        
        facet_wrap(~strain_year, ncol = 5) +
        
        coord_cartesian(ylim = c(0, 9)) +
        
        theme_bw()
      
      plot(p)
    }
    dev.off()
  }
  
}

