
suppressMessages(library(tidyverse))
suppressMessages(library(arrow))

suppressMessages(library(tidybayes))
suppressMessages(library(bayesplot))
suppressMessages(library(ggtext))

save_hdf5 <- function(data_list, filename) {
  if(file.exists(filename)) {
    file.remove(filename)
  }
  
  for (i in 1:length(data_list)) {
    rhdf5::h5write(data_list[[i]], filename, names(data_list)[[i]])
  }
  
  return(filename)
}


read_hdf5 <- function(filename) {
  item_names <- rhdf5::h5ls(filename)$name
  
  l <- map(item_names,
           function(item_name) {
             rhdf5::h5read(filename, item_name)
           }) %>%
    `names<-`(item_names)
  
  return(l)
}

read_model_data <- function(filename) {
  data <- read_hdf5(filename)
  
  types <- map_chr(data, ~ class(.)[1])
  
  for(i in 1:length(data)) {
    if(types[[i]] == "data.frame") {
      data[[i]] <- process_data_df(data[[i]], data$modelled_years)
    }
  }
  
  return(data)
}

process_data_df <- function(data_df, modelled_years) {
  col_names <- colnames(data_df)
  
  if("ix_t_obs" %in% col_names) {
    data_df <- data_df %>%
      mutate(year_observed = modelled_years[ix_t_obs])
  }
  
  if("ix_strain" %in% col_names) {
    data_df <- data_df %>%
      mutate(strain_year = modelled_years[ix_strain])
  }
  
  if("ix_t_birth" %in% col_names) {
    data_df <- data_df %>%
      mutate(ix_t_birth = if_else(ix_t_birth == 0, NA_integer_, ix_t_birth))
    
    if(!("year_of_birth" %in% col_names)) {
      data_df <- data_df %>% 
        mutate(year_of_birth = modelled_years[ix_t_birth])
    }
  }
  
  if("ix_t" %in% col_names) {
    data_df <- data_df %>%
      mutate(year = modelled_years[ix_t])
  }
  
  data_df <- data_df %>%
    as_tibble()
  
  return(data_df)
}

clean_chain <- function(chain_df) {
  chain_df %>%
    rename(.chain = chain, .iteration = iteration) %>% 
    mutate(#Chain = .chain, # Copy for whatever reason.
           .draw = (.iteration - min(.iteration)) + (.chain - 1) * (max(.iteration) - min(.iteration) + 1),
           .before = 3) %>%
    rename_with(function(x) str_remove(x, " ")) # Remove spaces from array indexing
}

read_fit_data <- function(filename, modelled_years) {
  data <- read_hdf5(filename)
  
  data$chain <- clean_chain(data$chain)
  data$ppd <- process_data_df(data$ppd, modelled_years)
  
  return(data)
}

read_chain <- function(filename) {
  clean_chain(read_parquet(filename))
}

render_quarto <- function(run_name, run_dir, model_data_file, chain_file, quarto_file) {
  file_out <- str_c("reports/", run_name, "_", today(), ".pdf")
  file_tmp <- str_c(run_name, ".pdf")
  
  suppressWarnings(file.remove(file_out))
  quarto::quarto_render(
    quarto_file,
    output_file = file_tmp,
    execute_params = list(model_data_file = model_data_file, chain_file = chain_file),
    quiet = TRUE
  )
  
  file.copy(file_tmp, file_out)
  file.remove(file_tmp)
}

render_quarto_sim_study <- function(run_name, chain_name) {
  render_quarto(
    str_c(run_name, "_", chain_name), 
    str_c("runs/", run_name, "/"),
    str_c("runs/", run_name, "/model_data.hdf5"),
    str_c("runs/", run_name, "/chain_", chain_name, ".parquet"),
    "R/sim_study_chain_report.qmd"
  )
}


render_quarto_data_study <- function(run_name, chain_name) {
  render_quarto(
    str_c(run_name, "_", chain_name), 
    str_c("runs/", run_name, "/"),
    str_c("runs/", run_name, "/model_data.hdf5"),
    str_c("runs/", run_name, "/chain_", chain_name, ".parquet"),
    "R/data_study_chain_report.qmd"
  )
}


get_strain_year_from_name <- function(strain_name) {
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


get_strain_year_from_name <- function(strain_name) {
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

make_kucharski_antigenic_distances <- function(modelled_years) {
  raw_strain_coords <- read_csv("input_data/kucharski_2018/datasets/antigenic_coords.csv") %>%
    rename(strain_name = viruses, Y = AG_x, X = AG_y) %>%
    mutate(strain_year = get_strain_year_from_name(strain_name)) %>%
    arrange(strain_year)
  
  fit_strain_coords <- generate_antigenic_map(raw_strain_coords, modelled_years)
  
  antigenic_distances <- fit_strain_coords %>%
    filter(strain_year %in% modelled_years) %>%
    arrange(strain_year) %>% 
    generate_antigenic_distances()
  
  return(antigenic_distances)
}


make_gam_antigenic_distances <- function(modelled_years) {
  
  strain_coords <- read_csv("input_data/kucharski_2018/datasets/antigenic_coords.csv") %>%
    rename(strain_name = viruses, y = AG_x, x = AG_y) %>%
    mutate(strain_year = get_strain_year_from_name(strain_name)) %>%
    arrange(strain_year) %>%
    mutate(t = strain_year - min(strain_year))
  
  library(mgcv)
  
  fit_x <- gam(x ~ s(t, k = 12), data = strain_coords)
  fit_y <- gam(y ~ s(t, k = 12), data = strain_coords)
  
  gratia::draw(fit_x, residuals = TRUE)
  gratia::draw(fit_y, residuals = TRUE)
  
  pred <- tibble(
    strain_year = modelled_years
  ) %>%
    mutate(t = strain_year - min(strain_year)) %>%
    
    mutate(x = predict(fit_x, newdata = .),
           y = predict(fit_y, newdata = .))
  
  antigenic_distances <- fit_strain_coords %>%
    filter(strain_year %in% modelled_years) %>%
    arrange(strain_year) %>% 
    generate_antigenic_distances()
  
  return(antigenic_distances)
}
