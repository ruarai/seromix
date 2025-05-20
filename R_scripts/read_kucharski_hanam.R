

# Note no difference from this to HaNam_data.Rdata I believe
load("input_data/HaNam_data_V1.RData")

modelled_years <- 1968:2012

process_individual_year <- function(data) {
  data %>% 
    t() %>%
    as_tibble() %>%
    `colnames<-`(c("year_observed", "observed_titre", "strain_year", "ix_sample"))
}

obs_df <- map(1:length(test.list), function(i) {
  map(test.list[[i]], process_individual_year) %>%
    bind_rows() %>% 
    mutate(ix_subject = i)
}) %>%
  bind_rows() %>%
  
  mutate(ix_subject = as.integer(ix_subject),
         ix_t_obs = match(year_observed, modelled_years),
         ix_strain = match(strain_year, modelled_years))  %>%
  
  select(ix_subject, ix_t_obs, ix_strain, year_observed, strain_year, observed_titre) %>%
  drop_na(observed_titre)

tar_source()
model_data_hanam <- read_model_data("runs/hanam_2018/model_data.hdf5")

obs_A <- model_data_hanam$observations %>%
  select(ix_subject, ix_strain, ix_t_obs, observed_titre) %>%
  mutate(across(everything(), c))


obs_B <- obs_df %>%
  select(ix_subject, ix_strain, ix_t_obs, observed_titre)

# If no differences, our data reading code is good
compare(obs_A, obs_B)




