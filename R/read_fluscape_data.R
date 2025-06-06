

read_fluscape_data_neuts <- function() {
  data_raw <- read_csv("input_data/kucharski_2018/datasets/Fluscape_SupplmentalDataS1.csv")
  
  clean_titre <- function(x) {
    # Note extra exp term here, taken from Kucharski code
    # Change of base to log2 titre
    round(log2(exp(as.numeric(x)) / 10) + 1, digits = 6)
    # equivalent to x log2(e) - log2(10) +1
  }
  
  modelled_years <- 1968:2009
  test_year <- 2009
  
  observations_df <- data_raw %>%
    rename(ix_subject = id,
           vaccination_status = is.vac,
           observed_titre = titers,
           strain_name = neut.against,
           strain_year = year.strain,
           ix_location = loc) %>%
    mutate(year_observed = test_year,
           observed_titre = clean_titre(observed_titre),
           ix_t_obs = match(year_observed, modelled_years),
           ix_strain = match(strain_year, modelled_years)) %>%
    
    mutate(across(c(ix_subject, ix_t_obs, ix_strain), as.integer)) %>%
    
    select(
      ix_subject, ix_t_obs, ix_strain,
      vaccination_status,
      age,
      year_observed,
      strain_name,
      strain_year,
      observed_titre
    ) %>%
    
    arrange(ix_subject, ix_t_obs, ix_strain)
  
  
  subject_birth_data <- observations_df %>%
    distinct(ix_subject, age) %>% 
    mutate(year_of_birth = test_year - age,
           ix_t_birth = match(year_of_birth, modelled_years),
           ix_t_birth = replace_na(ix_t_birth, 0)) %>%
    select(ix_subject, year_of_birth, ix_t_birth)
  
  antigenic_distances <- make_kucharski_antigenic_distances(modelled_years)
  
  model_data <- list(
    observations = observations_df,
    
    modelled_years = modelled_years,
    antigenic_distances = antigenic_distances,
    subject_birth_data = subject_birth_data
  )
  
  return(model_data)
}


read_fluscape_data_HI <- function() {
  data_raw <- read_csv("input_data/kucharski_2018/datasets/Fluscape_HI_data.csv")
  
  clean_titre <- function(x) {
    y <- round(log2(as.numeric(x) / 10) + 1, digits = 6)
    
    if_else(is.infinite(y), 0, y)
  }
  
  modelled_years <- 1968:2009
  test_year <- 2009
  
  observations_df <- data_raw %>%
    rename(age = Age) %>% 
    mutate(ix_subject = row_number()) %>% 
    pivot_longer(-c(age, ix_subject),
                 names_to = "strain_name",
                 values_to = "observed_titre") %>%
    mutate(year_observed = test_year,
           observed_titre = clean_titre(observed_titre),
           strain_year = as.numeric(str_extract(strain_name, "\\d{4}$")),
           ix_t_obs = match(year_observed, modelled_years),
           ix_strain = match(strain_year, modelled_years)) %>%
    
    mutate(across(c(ix_subject, ix_t_obs, ix_strain), as.integer)) %>%
    
    select(
      ix_subject, ix_t_obs, ix_strain,
      age,
      year_observed,
      strain_name,
      strain_year,
      observed_titre
    ) %>%
    
    arrange(ix_subject, ix_t_obs, ix_strain)
  
  subject_birth_data <- observations_df %>%
    distinct(ix_subject, age) %>% 
    mutate(year_of_birth = test_year - age,
           ix_t_birth = match(year_of_birth, modelled_years),
           ix_t_birth = replace_na(ix_t_birth, 0)) %>%
    select(ix_subject, year_of_birth, ix_t_birth)
  
  # antigenic_distances <- make_kucharski_antigenic_distances(modelled_years)
  
  # Manually copy across, so that they are identical
  antigenic_distances <- read_rds("input_data/kucharski_2018/antigenic_distances_kucharski.rds")
  
  model_data <- list(
    observations = observations_df,
    
    modelled_years = modelled_years,
    antigenic_distances = antigenic_distances,
    subject_birth_data = subject_birth_data
  )
  
  return(model_data)
}

