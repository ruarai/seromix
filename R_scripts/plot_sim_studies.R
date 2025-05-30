

library(targets)
library(tidyverse)

tar_source()

base_run_dir <- "runs/sim_study_hanam_2018_4"

get_ix_run <- function(x, base_run_dir) {
  str_extract(x, str_c("(?<=", base_run_dir, "/)\\d+(?=/)"))
}

get_chain_name <- function(x) {
  str_match(x, "chain_(.+)\\.parquet")[,2]
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
  mutate(ix_run = as.numeric(ix_run))

n_alive <- model_data_ls[[1]]$subject_birth_data %>%
  mutate(ix_t_birth = replace_na(ix_t_birth, 0)) %>%
  expand_grid(ix_t = 1:length(model_data_ls[[1]]$modelled_years)) %>%
  group_by(ix_t) %>%
  summarise(n_alive = sum(ix_t >= ix_t_birth))

simulation_metadata_sim_ar <- model_data_ls %>% 
  map(function(x) {
    x$infections %>%
      left_join(n_alive, by = join_by(ix_t)) %>%
      group_by(ix_t) %>% 
      summarise(p_inf = n() / n_alive[1])
  }) %>%
  bind_rows(.id = "ix_run") %>%
  mutate(ix_run = as.numeric(ix_run)) %>%
  left_join(simulation_metadata)


continuous_params <- model_data_ls[[1]]$continuous_params %>%
  pivot_longer(everything(), names_to = "variable", values_to = "true_value")



chain_files <- list.files(
  base_run_dir, pattern = "chain_*",
  recursive = TRUE, full.names = TRUE
)


summarise_chain <- function(chain_file, ix_run, chain_name) {
  parnames <- c("mu_long", "mu_short", "sigma_long", "sigma_short", "tau", "obs_sd", "omega")
  
  read_chain(chain_file) %>%
    mutate(ix_run = as.numeric(ix_run),
           chain_name = chain_name) %>%
    
    filter(.iteration > 40000) %>%
    
    select(.iteration, .chain, ix_run, chain_name, any_of(parnames)) %>%
    
    group_by(chain_name, ix_run, .chain) %>% 
    
    summarise_draws()
}

draw_summaries <- tibble(
  chain_file = chain_files,
  ix_run = get_ix_run(chain_files, base_run_dir),
  chain_name = get_chain_name(chain_files)
) %>%
  pmap(summarise_chain) %>%
  bind_rows()

plot_data_summaries <- simulation_metadata %>%
  right_join(draw_summaries)

plot_data_summaries %>%
  filter(endemic_mean_ar == 0.5,
         variable == "mu_long") %>% 
  ggplot() +
  
  geom_linerange(aes(xmin = q5, xmax = q95, y = factor(ix_run), colour = factor(.chain)),
                 position = position_dodge2(width = 0.2)) +
  geom_point(aes(x = mean, y = factor(ix_run), colour = factor(.chain)),
             position = position_dodge2(width = 0.2)) +
  
  geom_vline(aes(xintercept = true_value),
             linetype = "dashed",
             continuous_params %>% filter(variable == "mu_long")) +
  
  facet_wrap(~chain_name, ncol = 1, scales = "free_y") +
  
  theme(legend.position = "none")



plot_data_summaries %>%
  ggplot() +
  geom_linerange(aes(ymin = 1, ymax = ess_bulk, x = endemic_mean_ar),
                 position = position_dodge2(width = 0.5)) +
  
  facet_grid(cols = vars(chain_name), rows = vars(variable), scales = "free_y")


summarise_chain_attack_rate <- function(chain_file, ix_run, chain_name) {
  read_chain(chain_file) %>%
    
    filter(.iteration > 40000) %>%
    
    get_attack_rate(model_data_ls[[1]]$subject_birth_data) %>%
    
    mutate(ix_run = as.numeric(ix_run),
           chain_name = chain_name,
           .before = 1)
}

draw_summaries_ar <- tibble(
  chain_file = chain_files,
  ix_run = get_ix_run(chain_files, base_run_dir),
  chain_name = get_chain_name(chain_files)
) %>%
  pmap(summarise_chain_attack_rate) %>%
  bind_rows()


plot_data_summaries_ar <- simulation_metadata %>%
  right_join(draw_summaries_ar)


ggplot() +
  geom_line(aes(x = ix_t, y = p, group = interaction(ix_run, chain_name)),
            plot_data_summaries_ar) +
  
  geom_hline(aes(yintercept = endemic_mean_ar), tibble(endemic_mean_ar = unique(plot_data_summaries_ar$endemic_mean_ar)),
             linetype = "dashed") +
  
  facet_grid(rows = vars(endemic_mean_ar),
             cols = vars(chain_name))

