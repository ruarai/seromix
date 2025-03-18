
source("R/common.R")

hanam_data_raw <- read_csv("input_data/HaNamCohort.csv", col_types = cols(.default = col_character()))



n_strain <- ncol(hanam_data_raw) - 2

clean_hi_titer <- function(x) {
  case_when(
    x == "*" ~ NA_real_,
    TRUE ~ log2(as.numeric(x) / 10) + 1
  )
}

get_strain_year <- function(strain_name) {
  year_part <- strain_name %>%
    str_split("/") %>%
    map(last) %>%
    str_extract("^\\d{2,4}")
  
  case_when(
    nchar(year_part) == 2 & as.numeric(year_part) > 15 ~ str_c("19", year_part),
    nchar(year_part) == 2 & as.numeric(year_part) <= 15 ~ str_c("20", year_part),
    nchar(year_part) == 4 ~ year_part
  ) %>%
    as.numeric()
}



hanam_data_processed <- hanam_data_raw %>%
  rename(subject_id = `Subject number`,
         year_sampled = `Sample year`) %>%
  mutate(subject_id = as.integer(subject_id),
         year_sampled = as.integer(year_sampled)) %>%
  
  pivot_longer(cols = -c(subject_id, year_sampled),
               names_to = "strain_name",
               values_to = "titre") %>%
  
  mutate(titre = clean_hi_titer(titre),
         strain_year = get_strain_year(strain_name)) %>%
  
  drop_na(titre)


unique_sample_years <- sort(unique(hanam_data_processed$year_sampled))
unique_strain_years <- sort(unique(hanam_data_processed$strain_year))

sample_or_strain_years <- sort(unique(c(unique_sample_years, unique_strain_years))) %>%
  as.integer()

modelled_years <- min(sample_or_strain_years):max(sample_or_strain_years)

baseline_year <- modelled_years[1]


observations_df <- hanam_data_processed %>%
  
  mutate(ix_t_obs = match(year_sampled, modelled_years),
         ix_strain = match(strain_year, modelled_years)) %>%
  
  select(ix_subject = subject_id, ix_t_obs, ix_strain, observed_titre = titre) %>%
  
  mutate(across(c(ix_subject, ix_t_obs, ix_strain), as.integer))



antigenic_coords <- read_csv("input_data/antigenic_coords.csv") %>%
  rename(strain_name = viruses, Y = AG_x, X = AG_y) %>%
  mutate(strain_year = get_strain_year(strain_name)) %>%
  arrange(strain_year)



fit_coords <- generate_antigenic_map(antigenic_coords)

ggplot() +
  geom_point(aes(x = X, y = Y, colour = strain_year),
             antigenic_coords) +
  
  geom_path(aes(x = x, y = y, colour = strain_year),
            fit_coords) +
  
  geom_point(aes(x = x, y = y, colour = strain_year),
             fit_coords %>% filter(strain_year %in% modelled_years))


antigenic_distances <- fit_coords %>%
  filter(strain_year %in% modelled_years) %>%
  arrange(strain_year) %>% 
  generate_antigenic_distances()



subject_birth_data <- read_csv("input_data/HaNam_YOB.csv", col_names = "year_of_birth") %>%
  mutate(ix_subject = row_number()) %>%
  mutate(ix_t_birth = match(year_of_birth, modelled_years),
         ix_t_birth = replace_na(ix_t_birth, 0))


model_data <- list(
  observations_df = observations_df,
  modelled_years = modelled_years,
  antigenic_distances = antigenic_distances,
  subject_birth_ix = subject_birth_data$ix_t_birth
)

save_hdf5(model_data, "data/hanam_2018/model_data.hdf5")


write_csv(hanam_data_processed, "data/hanam_2018/data_processed.csv")











