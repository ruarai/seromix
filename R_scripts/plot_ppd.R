
library(targets)
library(tidyverse)

tar_source()

run_name <- "hanam_2018"
chain_name <- "linear_basic_2"

run_dir <- str_c("runs/", run_name, "/")
model_data <- read_model_data(str_c(run_dir, "/model_data.hdf5"))


ppd <- arrow::read_parquet(str_c(run_dir, "ppd_", chain_name, ".parquet")) |>
  process_data_df(model_data$modelled_years)


chain <- read_chain(str_c(run_dir, "chain_", chain_name, ".parquet"))

inf_data <- get_inf_data_data_study(chain |> filter(.iteration > 5000), model_data = model_data) |>
  process_data_df(model_data$modelled_years)

ix_subject_ex <- 5


# plot_data_ppd |>
#   ggplot() +
#   geom_tile(aes(x = strain_year, y = ix_draw, fill = observed_titre)) +
#   
#   facet_wrap(~year_observed) +
#   
#   scale_x_continuous(labels = function(x) { map_chr(x, function(y) str_c("'", str_sub(y, 3, 4)))}) +
#   
#   scale_fill_viridis_c()


subjects <- unique(ppd$ix_subject)

pdf(str_c(run_dir, "ppd_", chain_name, ".pdf"), width = 10, height = 8)
for(ix_subject_ex in subjects) {
  plot_data_ppd <- ppd |> 
    filter(ix_subject == ix_subject_ex,
           ix_draw <= 30)
  
  plot_data_obs <- model_data$observations |>
    filter(ix_subject == ix_subject_ex)
  
  plot_data_inf <- inf_data |> 
    mutate(year_observed = year) |>
    filter(ix_subject == ix_subject_ex)
  
  p <- ggplot() +
    
    geom_line(aes(x = strain_year, y = observed_titre, group = ix_draw),
              alpha = 0.5, linewidth = 0.5, colour = "grey30",
              plot_data_ppd) +
    
    geom_linerange(aes(x = year, ymin = 0, ymax = 8),
                   colour = "black", linewidth = 0.3,
                   plot_data_inf) +
    
    geom_linerange(aes(x = year, ymin = 0, ymax = 8 * p),
                   colour = "green", linewidth = 1.5,
                   plot_data_inf) +
    
    geom_point(aes(x = strain_year, y = observed_titre),
               colour = "red", size = 0.5,
               plot_data_obs) +
    
    scale_x_continuous(labels = function(x) { map_chr(x, function(y) str_c("'", str_sub(y, 3, 4)))}) +
    
    facet_wrap(~year_observed)
  
  plot(p)
}
dev.off()


# raw_strain_coords <- read_csv("input_data/kucharski_2018/datasets/antigenic_coords.csv") |>
#   rename(strain_name = viruses, Y = AG_x, X = AG_y) |>
#   mutate(strain_year = get_strain_year_from_name(strain_name)) |>
#   arrange(strain_year)
# 
# fit_strain_coords <- generate_antigenic_map(raw_strain_coords, model_data$modelled_years)
# 
# strain_locations_z <- fit_strain_coords |>
#   mutate(dist = sqrt((x-lag(x))^2 + (y - lag(y))^2),
#          dist = replace_na(dist, 0),
#          z = cumsum(dist))
# 
# plot_data_ppd <- ppd |> 
#   filter(ix_subject == ix_subject_ex,
#          ix_draw <= 30) |>
#   left_join(strain_locations_z)
# 
# plot_data_obs <- model_data$observations |>
#   filter(ix_subject == ix_subject_ex) |>
#   left_join(strain_locations_z)
# 
# p <- ggplot() +
#   
#   geom_line(aes(x = z, y = observed_titre, group = ix_draw),
#             alpha = 0.5, linewidth = 0.5, colour = "grey30",
#             plot_data_ppd) +
#   
#   geom_point(aes(x = z, y = observed_titre),
#              colour = "red", size = 0.5,
#              plot_data_obs) +
#   
#   facet_wrap(~year_observed)
# 
# plot(p)




# 
# ix_subject_ex <- 24
# 
# plot_data_ppd <- ppd |>
#   filter(ix_subject == ix_subject_ex,
#          ix_draw <= 30)
# 
# plot_data_obs <- model_data$observations |>
#   filter(ix_subject == ix_subject_ex)
# 
# 
# ggplot() +
#   geom_line(aes(x = year_observed, y = observed_titre, group = ix_draw),
#             alpha = 0.5, linewidth = 0.5,
#             plot_data_ppd) +
# 
#   geom_point(aes(x = year_observed, y = observed_titre),
#              colour = "red", size = 0.5,
#              plot_data_obs) +
# 
#   scale_x_continuous(labels = function(x) { map_chr(x, function(y) str_c("'", str_sub(y, 3, 4)))}) +
# 
#   coord_cartesian(xlim = c(2000, 2012)) +
# 
#   facet_wrap(~strain_year)
# 
# 
