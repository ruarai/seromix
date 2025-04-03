

load("Simulated_data_SIM_1.RData")

modelled_years <- 1968:2012

process_individual_year <- function(data) {
  data %>% 
    t() %>%
    as_tibble() %>%
    
    rename(year_observed = test.year,
           observed_titre = titredat,
           strain_year = strain_years,
           ix_sample = sample.index)

}

obs_df <- map(1:length(test.listSim), function(i) {
  map(test.listSim[[i]], process_individual_year) %>%
    bind_rows() %>% 
    mutate(ix_subject = i)
}) %>%
  bind_rows() %>%
  
  mutate(ix_subject = as.integer(ix_subject),
         ix_t_obs = match(year_observed, modelled_years),
         ix_strain = match(strain_year, modelled_years)) %>%
  
  select(ix_subject, ix_t_obs, ix_strain, year_observed, strain_year, observed_titre)


subject_birth_data <- tibble(ix_subject = 1:n_part, ix_t_birth = 0)

antigenic_distances <- as.matrix(dist(antigenic.map.in))


infections <- reshape2::melt(historytabSim, varnames = c("ix_subject", "ix_t"), value.name = "infected") %>%
  as_tibble() %>%
  mutate(year = modelled_years[ix_t]) %>%
  filter(infected == 1)

model_data <- list(
  observations = obs_df,
  
  modelled_years = modelled_years,
  antigenic_distances = antigenic_distances,
  subject_birth_data = subject_birth_data,
  
  infections = infections
)

tar_source()

save_hdf5(model_data, "runs/hanam_ak_1/model_data.hdf5")


sample_dates <- obs_df %>%
  distinct(ix_subject, year_observed)

strains_observed <- obs_df %>%
  distinct(ix_subject, strain_year)

ggplot() +

  geom_tile(aes(x = year, y = ix_subject),
            fill = "lightblue",
            infections) +
  
  geom_point(aes(x = year_observed, y = ix_subject),
             size = 0.5,
             sample_dates) +
  
  geom_rug(aes(x = year), tibble(year = modelled_years)) +
  
  theme_bw()

ggplot() +
  
  geom_tile(aes(x = year, y = ix_subject),
            fill = "lightblue",
            infections)  +
  
  geom_tile(aes(x = strain_year, y = ix_subject),
            fill = "red", alpha = 0.3,
            strains_observed) +
  geom_point(aes(x = year_observed, y = ix_subject),
             sample_dates) +
  
  geom_rug(aes(x = year), tibble(year = modelled_years)) +
  
  theme_bw()
