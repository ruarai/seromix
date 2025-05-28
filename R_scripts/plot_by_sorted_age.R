

model_data <- tar_read(hanam_data_age)

# age_data <- read_csv("input_data/kgostic_data_emporium/Fonville_2014/YOB_inferred.csv") %>%
#   select(ix_subject = `Subject number`,
#          year_of_birth = YOB) %>%
#   mutate(ix_t_birth = match(year_of_birth, modelled_years),
#          ix_t_birth = replace_na(ix_t_birth, 0),
#          ix_subject = as.integer(ix_subject)) %>%
#   mutate(ix_subject_sorted = 70 - row_number())


waldo::compare(
  age_data %>%
    distinct(ix_subject, year_of_birth) %>%
    arrange(ix_subject, year_of_birth) %>%
    mutate(across(everything(), as.integer)),
  age_data_2 %>%
    distinct(ix_subject, year_of_birth) %>%
    arrange(ix_subject, year_of_birth) %>%
    mutate(across(everything(), as.integer))
)

observations <- model_data$observations %>%
  left_join(age_data)

plot_dir <- "runs/hanam_2018_age/"

n_pages <- ceiling(max(observations$ix_subject) / 5)
pdf(str_c(plot_dir, "data_by_year_observed_sorted.pdf"), width = 10, height = 6)

year_blank <- tibble(year_observed = 2007:2012)


for(i in 1:n_pages) {
  p <- observations %>%
    filter(ix_subject_sorted %in% (((i - 1) * 5):(i * 5 - 1) + 1) ) %>% 
    ggplot() +
    geom_point(aes(x = strain_year, y = observed_titre),
               size = 0.7) +
    
    geom_blank(aes(y = 0), year_blank) +
    
    facet_grid(rows = vars(ix_subject_sorted), cols = vars(year_observed)) +
    
    coord_cartesian(ylim = c(0, 9)) +
    
    theme_bw()
  
  plot(p)
}
dev.off()


pdf(str_c(plot_dir, "data_by_year_observed.pdf"), width = 10, height = 6)
for(i in 1:n_pages) {
  p <- observations %>%
    filter(ix_subject %in% (((i - 1) * 5):(i * 5 - 1) + 1) ) %>% 
    ggplot() +
    geom_point(aes(x = strain_year, y = observed_titre),
               size = 0.7) +
    
    geom_blank(aes(y = 0), year_blank) +
    
    facet_grid(rows = vars(ix_subject), cols = vars(year_observed)) +
    
    coord_cartesian(ylim = c(0, 9)) +
    
    theme_bw()
  
  plot(p)
}
dev.off()
