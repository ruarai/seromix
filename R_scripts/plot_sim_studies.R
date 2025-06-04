

library(targets)
library(tidyverse)

tar_source()

base_run_dir <- "runs/sim_study_hanam_2018_4"


get_ix_run <- function(x, base_run_dir) {
  str_extract(x, str_c("(?<=", base_run_dir, "/)\\d+(?=/)"))
}

model_data_files <- list.files(
  base_run_dir, pattern = "model_data.hdf5",
  recursive = TRUE, full.names = TRUE
)
names(model_data_files) <- get_ix_run(model_data_files, base_run_dir)

model_data_ls <- model_data_files %>%
  map(read_model_data)


simulation_metadata <- model_data_ls %>%
  map(~ .x$simulation_metadata) %>%
  bind_rows(.id = "ix_run") %>%
  mutate(ix_run = as.numeric(ix_run),
         total_infections = map_dbl(model_data_ls, ~ nrow(.x$infections))) %>%
  
  filter(endemic_mean_ar %in% c(0.1, 0.3, 0.5))


total_time_points <- model_data_ls[[1]]$subject_birth_data %>%
  mutate(n_alive = 45 - replace_na(ix_t_birth, 0)) %>%
  summarise(n_alive = sum(n_alive)) %>%
  pull(n_alive)


draw_summaries <- list.files(
  base_run_dir, pattern = "^chain_.*summary\\.csv$*",
  recursive = TRUE, full.names = TRUE
) %>%
  map(read_csv, show_col_types = FALSE) %>%
  bind_rows() %>%
  rename(ix_run = run_ix)

continuous_params <- model_data_ls[[1]]$continuous_params %>%
  pivot_longer(everything(), names_to = "variable", values_to = "true_value")


plot_data_summaries <- simulation_metadata %>%
  right_join(draw_summaries) %>% 
  filter(endemic_mean_ar %in% c(0.1, 0.3, 0.5))



plot_data_subset <- plot_data_summaries %>% 
  filter(!str_detect(chain_name, "no_age"),
         !str_detect(chain_name, "uncorrected"))

variables <- unique(plot_data_summaries$variable)

x_var <- variables[1]

plot_data_subset %>%
  filter(variable == x_var) %>% 
  ggplot() +
  
  geom_linerange(aes(xmin = q95_lower, xmax = q95_upper, y = factor(ix_run %% 3), colour = factor(chain)),
                 position = position_dodge2(width = 0.2)) +
  geom_point(aes(x = mean, y = factor(ix_run %% 3), colour = factor(chain)),
             position = position_dodge2(width = 0.2)) +
  
  geom_vline(aes(xintercept = true_value),
             linetype = "dashed",
             continuous_params %>% filter(variable == x_var)) +
  
  facet_grid(rows = vars(chain_name), cols = vars(endemic_mean_ar),
             scales = "free_y", switch = "y") +
  
  xlab("Value") + ylab(NULL) +
  
  theme(legend.position = "none",
        strip.text.y.left = element_text(angle = 0),
        axis.text.y = element_blank()) +
  
  ggtitle(x_var)

plot_data_subset %>%
  filter(variable == "total_infections") %>% 
  ggplot() +
  
  geom_linerange(aes(xmin = q95_lower / total_time_points, xmax = q95_upper / total_time_points, y = factor(ix_run %% 3), colour = factor(chain)),
                 position = position_dodge2(width = 0.2)) +
  geom_point(aes(x = mean / total_time_points, y = factor(ix_run %% 3), colour = factor(chain)),
             position = position_dodge2(width = 0.2)) +
  
  
  geom_linerange(aes(x = total_infections / total_time_points, ymin = ix_run %% 3 + 0.5, ymax = ix_run %% 3 + 1.5),
                 simulation_metadata) +
  
  facet_grid(rows = vars(chain_name), cols = vars(endemic_mean_ar),
             switch = "y") +
  
  xlab("Value") + ylab(NULL) +
  
  theme(legend.position = "none",
        strip.text.y.left = element_text(angle = 0),
        axis.text.y = element_blank())
