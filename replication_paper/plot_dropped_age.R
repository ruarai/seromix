

library(targets)
tar_source()

source("replication_paper/common.R")

summary <- tar_read(combined_sim_ar_singular_summaries)

true_values <- tibble(
  mu_long = 2.0,
  mu_short = 2.0,
  omega = 0.75,
  sigma_long = 0.15,
  sigma_short = 0.05,
  tau = 0.05,
  obs_sd = 1.5
) |> 
  pivot_longer(everything(), names_to = "variable") |> 
  mutate(variable = factor(variable, var_names, var_labels))


prior_labels <- c(
  "BetaBernoulli_1_1" = "BetaBernoulli(1, 1)",
  "BetaBernoulliTimeVarying_1_1" = "BetaBernoulliTV(1, 1)",
  "BetaBernoulliSubjectVarying_1_1" = "BetaBernoulliSV(1, 1)"
)


plot_data <- summary |> 
  
  filter(variable %in% var_names,
         variable != "total_inf") |> 
  mutate(variable = factor(variable, var_names, var_labels),
         prior_name = prior_labels[prior_description],
         endemic_mean_ar = factor(endemic_mean_ar))


plot_data |> 
  # filter(prior_description == "BetaBernoulliTimeVarying_1_1") |>
  filter(prior_description == "BetaBernoulliSubjectVarying_1_1") |>
  # filter(prior_description == "BetaBernoulli_1_1") |> 
  ggplot() +
  
  geom_vline(aes(xintercept = value),
             linetype = "14",
             true_values) +
  
  geom_point(aes(x = median, y = endemic_mean_ar, colour = drop_age),
             size = 0.8,
             position = position_dodge2(width = 0.3)) +
  
  geom_linerange(aes(xmin = q95_lower, xmax = q95_upper, y = endemic_mean_ar, colour = drop_age),
                 position = position_dodge2(width = 0.3)) +
  
  facet_wrap(~variable, ncol = 4, scales = "free_x") +
  
  ggokabeito::scale_colour_okabe_ito(order = c(9, 5), name = "Age", labels = c("Kept", "Dropped")) +
  
  xlab("Value (median, 95% CrI)") + ylab(NULL) +
  plot_theme_paper +
  
  theme(strip.text = element_markdown(size = 12, family = "Utopia"),
        panel.grid.major.y = element_gridline,
        panel.spacing.x = unit(1, "cm"),
        legend.position = "bottom")


plot_data |> 
  filter(variable %in% var_labels[c(1, 3, 5)]) |> 
  ggplot() +
  
  geom_vline(aes(xintercept = value),
             linetype = "14",
             true_values |> filter(variable  %in% var_labels[c(1, 3, 5)])) +
  
  geom_point(aes(x = median, y = endemic_mean_ar, colour = drop_age),
             size = 0.8,
             position = position_dodge2(width = 0.3)) +
  
  geom_linerange(aes(xmin = q95_lower, xmax = q95_upper, y = endemic_mean_ar, colour = drop_age),
                 position = position_dodge2(width = 0.3)) +
  
  facet_grid(cols = vars(variable), 
             rows = vars(prior_name), scales = "free_x") +
  
  ggokabeito::scale_colour_okabe_ito(order = c(9, 5), name = "Age", labels = c("Kept", "Dropped")) +
  
  xlab("Value (median, 95% CrI)") + ylab(NULL) +
  plot_theme_paper +
  
  theme(strip.text = element_markdown(size = 12, family = "Utopia"),
        panel.grid.major.y = element_gridline,
        panel.background = element_rect(fill = "transparent", colour = "grey80"),
        panel.spacing.x = unit(1, "cm"),
        legend.position = "bottom")
