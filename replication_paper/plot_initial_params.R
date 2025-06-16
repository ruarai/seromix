
library(targets)
tar_source()

source("replication_paper/common.R")

summary <- tar_read(combined_summaries) |>
  filter(exp_group == "initial_conditions")

chains <- tar_read(combined_chains) |>
  filter(exp_group == "initial_conditions")

name_order <- c(
  "kucharski_data_study",
  "broad"
)

name_labels <- c(
  "Narrow",
  "Broad"
)

plot_data <- summary |>
  mutate(name = initial_params_name) |> 
  filter(run_name == "hanam_2018",
         prior_description == "BetaBernoulli_1_1",
         variable %in% var_names) |>
  mutate(name = fct_rev(factor(name, name_order, name_labels)),
         variable = factor(variable, var_names, var_labels))

ggplot(plot_data) +
  
  geom_point(aes(x = median, y = name),
             size = 0.8,
             position = position_dodge2(width = 0.3)) +
  geom_linerange(aes(xmin = q95_lower, xmax = q95_upper, y = name),
                 position = position_dodge2(width = 0.3)) +
  
  geom_vline(aes(xintercept = median),
             linetype = "14",
             plot_data |> filter(name == "Kucharski 2018")) +
  
  facet_wrap(~variable,
             ncol = 4, scales = "free_x") +
  
  xlab("Value (median, 95% CrI)") + ylab(NULL) +
  
  plot_theme_paper +
  
  theme(strip.text = element_markdown(size = 12, family = "Utopia"),
        panel.grid.major.y = element_gridline)

plot_data <- chains |>
  mutate(name = initial_params_name) |> 
  filter(run_name == "hanam_2018",
         prior_description == "BetaBernoulli_1_1") |>
  mutate(name = fct_rev(factor(name, name_order, name_labels))) |> 
  filter(.iteration > 150000)


ggplot(plot_data) +
  geom_point(aes(x = tau, y = mu_long, colour = factor(.chain))) +
  
  facet_wrap(~name) +
  
  plot_theme_paper +
  
  theme(panel.grid.major = element_gridline)


chain <- tar_read(chain_hanam_2018_11)

chain_inf <- chain |>
  clean_chain() |> 
  add_total_infections(tar_read(hanam_2018))

chain_inf |> 
  filter(.iteration > 50000) |>
  ggplot() +
  geom_point(aes(x = total_inf, y = obs_sd, colour = factor(.chain)),
             size = 0.5) +
  plot_theme_paper +
  
  theme(panel.grid.major = element_gridline)


chain_inf |> 
  filter(.iteration > 50000) |>
  ggplot() +
  geom_point(aes(x = total_inf, y = mu_long, colour = factor(.chain)),
             size = 0.5) +
  plot_theme_paper +
  
  theme(panel.grid.major = element_gridline)


model_data <- tar_read(hanam_2018)

col_names <- colnames(chain)
inf_names <- col_names[str_starts(col_names, "infections")]


results <- chain |> 
  clean_chain() |> 
  filter(.iteration > 50000,
         # .chain %in% c(6, 3, 2)
         ) |> 
  group_by(.chain) |> 
  summarise(across(starts_with("infections"), mean)) |>
  pivot_longer(starts_with("infections")) |>
  
  mutate(index = str_remove(name, "infections\\["),
         index = str_remove(index, "\\]"),
         index = str_split_fixed(index, ",", 2),
         
         ix_t = as.numeric(index[,1]), ix_subject = as.numeric(index[,2]))


ggplot(results) +
  geom_tile(aes(x = ix_t, y = ix_subject, fill = value)) +
  
  scale_fill_viridis_c() +
  
  facet_wrap(~.chain)


results |> 
  group_by(.chain, ix_t) |> 
  summarise(value = sum(value)) |>
  mutate(value_cumulative = cumsum(value)) |> 

  ggplot() +
  geom_step(aes(x = ix_t, y = value_cumulative, colour = factor(.chain))) +
  
  plot_theme_paper


library(loo)

loo::extract_log_lik(chain)




