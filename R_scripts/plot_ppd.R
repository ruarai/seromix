
tar_source()

run_name <- "hanam_2018"
chain_name <- "linear_basic_2"

run_dir <- str_c("runs/", run_name, "/")
model_data <- read_model_data(str_c(run_dir, "/model_data.hdf5"))


ppd <- arrow::read_parquet(str_c(run_dir, "ppd_", chain_name, ".parquet")) %>%
  process_data_df(model_data$modelled_years)


ix_subject_ex <- 16


subjects <- unique(ppd$ix_subject)

pdf(str_c(run_dir, "ppd_", chain_name, ".pdf"), width = 10, height = 8)
for(ix_subject_ex in subjects) {
  plot_data_ppd <- ppd %>% 
    filter(ix_subject == ix_subject_ex,
           ix_draw <= 30)
  
  plot_data_obs <- model_data$observations %>%
    filter(ix_subject == ix_subject_ex)
  
  
  p <- ggplot() +
    geom_line(aes(x = strain_year, y = observed_titre, group = ix_draw),
              alpha = 0.5, linewidth = 0.5,
              plot_data_ppd) +
    
    geom_point(aes(x = strain_year, y = observed_titre),
               colour = "red", size = 0.5,
               plot_data_obs) +
    
    scale_x_continuous(labels = function(x) { map_chr(x, function(y) str_c("'", str_sub(y, 3, 4)))}) +
    
    facet_wrap(~year_observed)
  
  plot(p)
}
dev.off()

obs_years <- unique(model_data$observations$year_observed)
pdf(str_c(run_dir, "ppd_error_", chain_name, ".pdf"), width = 10, height = 8)
for(ix_subject_ex in subjects) {
  
  plot_data_obs <- model_data$observations %>%
    filter(ix_subject == ix_subject_ex) %>%
    rename(observed_titre_true = observed_titre)
  
  plot_data_ppd <- ppd %>% 
    filter(ix_subject == ix_subject_ex,
           ix_draw <= 30) %>% 
    inner_join(plot_data_obs,
               by = join_by(ix_t_obs, ix_strain, ix_subject, year_observed, strain_year))
  
  p <- ggplot() +
    
    geom_point(aes(x = strain_year, y = observed_titre - observed_titre_true),
               size = 0.5, position = position_jitter(width = 0.5), alpha = 0.5,
               plot_data_ppd) +
    
    geom_blank(aes(x = 0), tibble(year_observed = obs_years)) +
    
    geom_hline(yintercept = 0) +
    
    scale_x_continuous(labels = function(x) { map_chr(x, function(y) str_c("'", str_sub(y, 3, 4)))}) +
    
    facet_wrap(~year_observed)
  
  plot(p)
}
dev.off()


