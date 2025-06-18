
library(targets)
tar_source()

source("replication_paper/common.R")

summary <- tar_read(combined_summaries) |>
  filter(exp_group == "initial_conditions")

chains <- tar_read(combined_chains) |>
  filter(exp_group == "initial_conditions")

var_labels_asterisk <- c(
  "<i>μ</i><sub>L</sub> (long-term boost)", 
  "<i>μ</i><sub>S</sub> (short-term boost)",
  "<i>σ</i><sub>L</sub> (long-term cross<br>reactivity) *", 
  "<i>σ</i><sub>S</sub> (short-term cross<br>reactivity) *",
  "<i>τ</i> (seniority) *",
  "<i>ω</i> (waning)",
  "<i>σ</i><sub>obs</sub> (obs sd.)", 
  "Total infections"
)

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
         variable_label = factor(variable, var_names, var_labels),
         variable_asterisk = factor(variable, var_names, var_labels_asterisk))

ggplot(plot_data) +
  
  geom_point(aes(x = median, y = name),
             size = 0.8,
             position = position_dodge2(width = 0.3)) +
  geom_linerange(aes(xmin = q95_lower, xmax = q95_upper, y = name),
                 position = position_dodge2(width = 0.3)) +
  
  facet_wrap(~variable_label,
             ncol = 4, scales = "free_x") +
  
  xlab("Value (median, 95% CrI)") + ylab(NULL) +
  
  plot_theme_paper +
  
  theme(strip.text = element_markdown(size = 12, family = "Utopia"),
        panel.grid.major.y = element_gridline)


ggsave(
  "replication_paper/results/initial_params.png",
  device = png,
  width = 12,
  height = 4, bg = "white"
)


ggplot(plot_data) +
  
  geom_point(aes(x = median, y = name),
             size = 0.8,
             position = position_dodge2(width = 0.3)) +
  geom_linerange(aes(xmin = q95_lower, xmax = q95_upper, y = name),
                 position = position_dodge2(width = 0.3)) +
  
  facet_wrap(~variable_asterisk,
             ncol = 4, scales = "free_x") +
  
  xlab("Value (median, 95% CrI)") + ylab(NULL) +
  
  plot_theme_paper +
  
  ggh4x::scale_x_facet(
    variable_asterisk == var_labels_asterisk[[5]],
    limits = c(0, 0.1)
  ) +
  
  ggh4x::scale_x_facet(
    variable_asterisk == var_labels_asterisk[[3]],
    limits = c(0.07, 0.2)
  ) +
  
  ggh4x::scale_x_facet(
    variable_asterisk == var_labels_asterisk[[4]],
    limits = c(0, 0.1)
  ) +

  theme(strip.text = element_markdown(size = 12, family = "Utopia"),
      panel.grid.major.y = element_gridline)


ggsave(
  "replication_paper/results/initial_params_zoomed.png",
  device = png,
  width = 12,
  height = 4, bg = "white"
)



plot_data <- chains |>
  mutate(name = initial_params_name) |> 
  # filter(proposal_name == "kucharski_data_study") |>
  mutate(name = fct_rev(factor(name, name_order, name_labels))) |> 
  filter(name == "Narrow") |> 
  filter(.iteration > 150000)


ggplot(plot_data) +
  geom_point(aes(x = tau, y = total_inf, colour = factor(.chain))) +
  
  facet_wrap(~name) +
  
  plot_theme_paper +
  
  coord_cartesian(xlim = c(0, 0.1),
                  ylim = c(0, 800)) +
  
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




