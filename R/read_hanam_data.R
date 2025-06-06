




read_hanam_data <- function(use_inferred_age = FALSE) {
  
  hanam_data_raw <- read_csv("input_data/kucharski_2018/datasets/HaNamCohort.csv", col_types = cols(.default = col_character()))
  
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
  
  
  if(use_inferred_age) {
    subject_birth_data <- read_csv("input_data/kgostic_data_emporium/Fonville_2014/YOB_inferred.csv") %>%
      select(ix_subject = `Subject number`,
             year_of_birth = YOB) %>%
      mutate(ix_t_birth = match(year_of_birth, modelled_years),
             ix_t_birth = replace_na(ix_t_birth, 0),
             ix_subject = as.integer(ix_subject)) %>%
      arrange(ix_subject)
  } else {
    n_subjects <- length(unique(observations_df$ix_subject))
    # As in Kucharski (2018), assume no known ages
    subject_birth_data <- tibble(
      ix_subject = as.integer(1:n_subjects),
      ix_t_birth = 0,
      year_of_birth = modelled_years[1] - 1 # Necessary for plots
    )
  }
  
  # Take age distribution
  # Likely from Fonville (2014)
  age_distribution <- read_csv(
    "input_data/kucharski_2018/datasets/HaNam_YOB.csv",
    col_names = "age", col_types = cols(age = col_integer())
  )$age

  
  # antigenic_distances <- make_kucharski_antigenic_distances(modelled_years)
  
  # Manually copy across, so that they are identical
  antigenic_distances <- read_rds("input_data/kucharski_2018/antigenic_distances_kucharski.rds")
  
  model_data <- list(
    observations = observations_df,
    age_distribution = age_distribution,
    
    modelled_years = modelled_years,
    antigenic_distances = antigenic_distances,
    subject_birth_data = subject_birth_data,
    
    initial_infections_manual = initial_infections
  )
  
  
  return(model_data)
}
