




read_hanam_data <- function() {
  
  hanam_data_raw <- read_csv("input_data/kucharski_2018/datasets/HaNamCohort.csv", col_types = cols(.default = col_character()))
  birth_data_raw <- read_csv("input_data/kucharski_2018/datasets/HaNam_YOB.csv", col_names = "year_of_birth")
  
  # Per Kucharski model, specify initial infections matrix
  initial_infections <- read_csv("input_data/kucharski_2018/R_datasets/hist_IC_H3.csv") %>%
    select(-1) %>%
    as.matrix() %>%
    t()
  
  
  
  
  n_strain <- ncol(hanam_data_raw) - 2
  
  clean_hi_titer <- function(x) {
    case_when(
      x == "*" ~ NA_real_,
      TRUE ~ suppressWarnings(log2(as.numeric(x) / 10) + 1)
    )
  }
  
  hanam_data_processed <- hanam_data_raw %>%
    rename(ix_subject = `Subject number`,
           year_observed = `Sample year`) %>%
    mutate(ix_subject = as.integer(ix_subject),
           year_observed = as.integer(year_observed)) %>%
    
    pivot_longer(cols = -c(ix_subject, year_observed),
                 names_to = "strain_name",
                 values_to = "observed_titre") %>%
    
    mutate(observed_titre = clean_hi_titer(observed_titre),
           strain_year = get_strain_year_from_name(strain_name)) %>%
    
    drop_na(observed_titre)
  
  
  unique_sample_years <- sort(unique(hanam_data_processed$year_observed))
  unique_strain_years <- sort(unique(hanam_data_processed$strain_year))
  
  sample_or_strain_years <- sort(unique(c(unique_sample_years, unique_strain_years))) %>%
    as.integer()
  
  modelled_years <- min(sample_or_strain_years):max(sample_or_strain_years)
  
  baseline_year <- modelled_years[1]
  
  
  observations_df <- hanam_data_processed %>%
    
    mutate(ix_t_obs = match(year_observed, modelled_years),
           ix_strain = match(strain_year, modelled_years)) %>%
    
    mutate(across(c(ix_subject, ix_t_obs, ix_strain), as.integer))
  
  
  # No age data available for this dataset
  subject_birth_data <- birth_data_raw %>%
    mutate(ix_subject = row_number()) %>% 
    mutate(ix_t_birth = 0)
  
  antigenic_distances <- make_kucharski_antigenic_distances(modelled_years)
  
  
  model_data <- list(
    observations = observations_df,
    
    modelled_years = modelled_years,
    antigenic_distances = antigenic_distances,
    subject_birth_data = subject_birth_data,
    
    initial_infections_manual = initial_infections
  )
  
  return(model_data)
}
