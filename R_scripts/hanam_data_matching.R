
model_data <- tar_read(hanam_data_age)

observations <- model_data$observations %>%
  left_join(age_data)


observations %>%
  group_by(ix_subject) %>%
  summarise(n_years_observed = length(unique(year_observed)))


obs_table <- observations %>%
  group_by(ix_subject, year_observed) %>%
  summarise(n_strains_observed = length(unique(strain_name))) %>%
  
  pivot_wider(names_from = year_observed, values_from = n_strains_observed,
              names_prefix = "obs_")


obs_table <- observations %>%
  group_by(ix_subject, year_observed) %>%
  summarise(extent_observed = max(strain_year) - min(strain_year)) %>%

  mutate(extent_observed = case_when(
    extent_observed == 43 ~ "full",
    TRUE ~ "partial",
  )) %>%
  
  pivot_wider(names_from = year_observed, values_from = extent_observed,
              names_prefix = "obs_")

obs_table %>%
  ungroup() %>% 
  count(obs_2007, obs_2008, obs_2009, obs_2010, obs_2011, obs_2012) %>%
  View()


fonville_age_table <- read_csv("input_data/fonville_age_matching.csv")


fonville_age_table %>%
  # select(-obs_2012) %>% 
  left_join(obs_table, relationship = "many-to-many") %>%
  
  select(ix_sorted, ix_subject) %>%
  group_by(ix_sorted) %>%
  summarise(ix_subject = str_c(ix_subject, collapse = ", ")) %>%
  
  View()


filter_matches <- function(tbl, row_candidate) {
  cols_to_check <- c("obs_2007", "obs_2008", "obs_2009", "obs_2010", "obs_2011")
  
  for (col_name in cols_to_check) {
    if (is.na(row_candidate[[col_name]])) {
      tbl <- tbl %>%
        filter(is.na(.data[[col_name]]))
    } else {
      tbl <- tbl %>%
        filter(!is.na(.data[[col_name]]))
    }
  }
  return(tbl)
}

ix_row <- 1

fonville_age_table[ix_row,]

possible_matches <- list()

for(ix_row in 1:44) {
  ix_matches <- obs_table %>% 
    filter_matches(fonville_age_table[ix_row,]) %>%
    pull(ix_subject)
  
  possible_matches <- c(possible_matches, list(ix_matches))
}


