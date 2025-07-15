

library(targets)
library(patchwork)
tar_source()

source("replication_paper/common.R")

summary_sim <- tar_read(combined_sim_ar_summaries)

summary_data <- tar_read(combined_summaries) |>
  filter(exp_group == "age_inclusion")

sim_true_values <- tibble(
  mu_long = 2.0,
  mu_short = 2.0,
  omega = 0.75,
  sigma_long = 0.15,
  sigma_short = 0.05,
  tau = 0.05,
  obs_sd = 1.5
) |> 
  pivot_longer(everything(), names_to = "variable") |> 
  mutate(variable = factor(variable, names(var_labels), var_labels))


prior_labels <- c(
  "BetaBernoulli_1_1" = "BetaBernoulli(1, 1)",
  "BetaBernoulliTimeVarying_1_1" = "BetaBernoulliTV(1, 1)",
  "BetaBernoulliSubjectVarying_1_1" = "BetaBernoulliSV(1, 1)"
)


sim_plot_data <- summary_sim |> 
  
  filter(variable %in% names(var_labels)[1:6],
         variable != "total_inf") |> 
  mutate(variable = factor(variable, names(var_labels), var_labels),
         prior_name = prior_labels[prior_description],
         endemic_mean_ar = factor(endemic_mean_ar))



p_sim <- sim_plot_data |> 
  filter(prior_name == "BetaBernoulli(1, 1)") |> 
  ggplot() +
  
  geom_vline(aes(xintercept = value),
             linetype = "14",
             sim_true_values |> filter(variable  %in% var_labels[1:6])) +
  
  geom_point(aes(x = median, y = endemic_mean_ar, colour = drop_age),
             size = 1.2,
             position = position_dodge2(width = 0.3)) +
  
  geom_linerange(aes(xmin = q95_lower, xmax = q95_upper, y = endemic_mean_ar, colour = drop_age),
                 linewidth = 0.7,
                 position = position_dodge2(width = 0.3)) +
  
  facet_wrap(~variable, scales = "free_x") +
  
  ggokabeito::scale_colour_okabe_ito(order = c(9, 5), name = "Age", labels = c("Kept", "Dropped")) +
  
  xlab("Value (median, 95% CrI)") + ylab(NULL) +
  plot_theme_paper +
  
  theme(strip.text = element_markdown(size = 12, family = "Utopia"),
        panel.grid.major.y = element_gridline,
        panel.background = element_facet_background,
        panel.spacing.x = unit(1, "cm"),
        legend.position = "none") +
  
  ggtitle("Simulation study")


data_run_labels <- c(
  "age_inclusion_hanam_2018_1" = "Dropped",
  "age_inclusion_hanam_2018_age_1" = "Kept"
)


plot_data <- summary_data |> 
  
  filter(variable %in% names(var_labels)[1:6],
         variable != "total_inf") |> 
  
  mutate(variable = factor(variable, names(var_labels), var_labels),
         name = factor(name, names(data_run_labels), data_run_labels))


p_data <- plot_data |> 
  ggplot() +
  
  geom_point(aes(x = median, y = name, colour = name),
             size = 1.2,
             position = position_dodge2(width = 0.3)) +
  
  geom_linerange(aes(xmin = q95_lower, xmax = q95_upper, y = name, colour = name),
                 linewidth = 0.7,
                 position = position_dodge2(width = 0.3)) +
  
  facet_wrap(~variable, scales = "free_x") +
  
  ggokabeito::scale_colour_okabe_ito(order = c(9, 5), name = "Age", labels = c("Kept", "Dropped")) +
  
  xlab("Value (median, 95% CrI)") + ylab(NULL) +
  plot_theme_paper +
  
  theme(strip.text = element_markdown(size = 12, family = "Utopia"),
        panel.grid.major.y = element_gridline,
        panel.background = element_facet_background,
        panel.spacing.x = unit(1, "cm"),
        legend.position = "bottom") +
  
  ggtitle("Ha Nam study")

p_sim / p_data


ggsave(
  "replication_paper/results/age_inclusion.png",
  device = png,
  width = 12,
  height = 8, bg = "white"
)
