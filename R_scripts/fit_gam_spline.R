


tar_source()

strain_coords <- read_csv("input_data/kucharski_2018/datasets/antigenic_coords.csv") |>
  rename(strain_name = viruses, y = AG_x, x = AG_y) |>
  mutate(strain_year = get_strain_year_from_name(strain_name)) |>
  arrange(strain_year) |>
  mutate(t = strain_year - min(strain_year))

library(mgcv)

fit_x <- gam(x ~ s(t, k = 12), data = strain_coords)
fit_y <- gam(y ~ s(t, k = 12), data = strain_coords)

gratia::draw(fit_x, residuals = TRUE)
gratia::draw(fit_y, residuals = TRUE)

pred <- tibble(
  strain_year = 1968:2012
) |>
  mutate(t = strain_year - min(strain_year)) |>
  
  mutate(x = predict(fit_x, newdata = .),
         y = predict(fit_y, newdata = .))


antigenic_distances <- pred |>
  arrange(strain_year) |> 
  generate_antigenic_distances()


antigenic_distances <- fit_strain_coords |>
  arrange(strain_year) |> 
  generate_antigenic_distances()

ggplot() +
  geom_point(aes(x = x, y = y, colour = strain_year),
             strain_coords) +
  geom_path(aes(x = x, y = y, colour = strain_year),
            pred) +
  geom_point(aes(x = x, y = y, colour = strain_year),
            pred)

ggplot() +
  geom_path(aes(x = x, y = y),
            pred) +
  geom_point(aes(x = x, y = y),
             pred) +
  
  geom_line(aes(x = x, y = y),
            colour = "grey50",
            fit_strain_coords) +
  
  geom_point(aes(x = x, y = y),
             colour = "grey50",
            fit_strain_coords)
