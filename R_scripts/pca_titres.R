


model_data <- tar_read(hanam_data)

observations <- model_data$observations


strain_names <- unique(observations$strain_name)
strain_names_years <- get_strain_year_from_name(strain_names)
strain_ordering <- order(strain_names_years)

strain_names[strain_ordering]
ordered_strain_names <- strain_names[strain_ordering]

observations <- model_data$observations %>%
  mutate(strain_name = factor(strain_name, ordered_strain_names))

observations %>%
  filter(ix_subject == 37) %>%
  
  ggplot() +
  geom_point(aes(x = year_observed, y = observed_titre)) +
  
  facet_wrap(~strain_year, ncol = 4)


observations %>%
  filter(ix_subject == 37) %>%
  
  ggplot() +
  geom_point(aes(x = year_observed, y = observed_titre, colour = factor(strain_year %% 4))) +
  
  facet_wrap(~strain_name, ncol = 4) +
  
  theme(legend.position = "none")



df_wide <- observations %>%
  select(ix_subject, year_observed, strain_name, observed_titre) %>%
  arrange(strain_name) %>% 
  pivot_wider(names_from = "strain_name", values_from = "observed_titre")


mat <- df_wide %>% 
  select(-c(ix_subject, year_observed)) %>% 
  as.matrix()

mat %>%
  exp() %>% 
  cor(use = "pairwise.complete.obs") %>%
  as_tibble(rownames = "strain_y") %>%
  pivot_longer(-c(strain_y), names_to = "strain_x") %>%
  
  mutate(strain_x = factor(strain_x, ordered_strain_names),
         strain_y = factor(strain_y, ordered_strain_names)) %>% 
  
  # group_by(strain_x) %>%
  # filter(!any(is.na(value))) %>% 
  
  ggplot() +
  geom_tile(aes(x = strain_x, y = strain_y, fill = value)) +
  
  scale_fill_viridis_c()


library(ade4)
pca <- nipals(2 ^ mat, tol = 1e-5, niter = 1000)

pca$co %>%
  as_tibble(rownames = "strain_name") %>%
  ggplot() +
  geom_point(aes(x = V1, y = V2, name = strain_name))

plotly::ggplotly()

dist_pca_long <- as.matrix(dist(pca$co)) %>%
  `rownames<-`(ordered_strain_names) %>%
  `colnames<-`(ordered_strain_names) %>% 
  as_tibble(rownames = "strain_y") %>%
  pivot_longer(-c(strain_y), names_to = "strain_x") %>%
  
  mutate(strain_x = factor(strain_x, ordered_strain_names),
         strain_y = factor(strain_y, ordered_strain_names))

dist_pca_long %>%
  
  ggplot() +
  geom_tile(aes(x = strain_x, y = strain_y, fill = value)) +
  
  scale_fill_viridis_c()


coords <- read_csv("input_data/kucharski_2018/datasets/antigenic_coords.csv")

coords_mat <- coords %>%
  rename(strain_name = viruses, x = AG_x, y = AG_y) %>%
  mutate(strain_name = factor(strain_name, ordered_strain_names)) %>%
  drop_na(strain_name) %>%
  select(x, y)

dist_coords_long <- as.matrix(dist(coords_mat)) %>%
  `rownames<-`(ordered_strain_names) %>%
  `colnames<-`(ordered_strain_names) %>% 
  as_tibble(rownames = "strain_y") %>%
  pivot_longer(-c(strain_y), names_to = "strain_x") %>%
  
  mutate(strain_x = factor(strain_x, ordered_strain_names),
         strain_y = factor(strain_y, ordered_strain_names))

dist_coords_long %>%
  
  ggplot() +
  geom_tile(aes(x = strain_x, y = strain_y, fill = value)) +
  
  scale_fill_viridis_c()


dist_coords_long %>%
  rename(dist_coords = value) %>% 
  left_join(dist_pca_long %>% rename(dist_pca = value))  %>%
  
  ggplot() +
  geom_point(aes(x = dist_coords, y = dist_pca))

