

e <- new.env()
load("input_data/kucharski_2018/R_datasets/FluScapeH3_HI_data.RData", envir = e)

modelled_years <- 1968:2009

process_individual_year <- function(data) {
  data |> 
    t() |>
    as_tibble() |>
    `colnames<-`(c("year_observed", "observed_titre", "strain_year", "ix_sample", "age"))
}

obs_df <- map(1:length(e$test.list), function(i) {
  map(e$test.list[[i]], process_individual_year) |>
    bind_rows() |> 
    mutate(ix_subject = i)
}) |>
  bind_rows() |>
  
  mutate(ix_subject = as.integer(ix_subject),
         ix_t_obs = match(year_observed, modelled_years),
         ix_strain = match(strain_year, modelled_years))  |>
  
  select(ix_subject, ix_t_obs, ix_strain, age, year_observed, strain_year, observed_titre) |>
  drop_na(observed_titre)

tar_source()



model_data <- read_fluscape_data_HI()

obs_A <- model_data$observations |>
  select(ix_subject, ix_strain, ix_t_obs, age, observed_titre) |>
  mutate(across(everything(), c)) |>
  arrange(ix_subject, ix_strain)


obs_B <- obs_df |>
  select(ix_subject, ix_strain, ix_t_obs, age, observed_titre) |>
  arrange(ix_subject, ix_strain)

# If no differences, our data reading code is good
compare(obs_A, obs_B)




