
library(targets)
tar_source()

source("replication_paper/common.R")

summary <- tar_read(combined_summaries) |>
  filter(exp_group == "model_comparison", mixture_importance_sampling)

lp_mixis_summary <- tar_read(combined_lp_mixis) %>%
  drop_na(lp_mixis) %>%
  mutate(variable = "mixis_lp", median = lp_mixis)

exp_labels <- c(
  "age_effect_2" = "Seniority and age-effect",
  "age_effect" = "Age-effect, no seniority",
  "kucharski" = "Seniority, no age-effect",
  "no_tau" = "No seniority or age-effect",
  "intercept" = "Seniority w/ intercept"
)

plot_data <- summary |>
  
  bind_rows(lp_mixis_summary) |> 
  mutate(name = exp_name) |>
  mutate(variable = factor(variable, names(var_labels), var_labels),
         name = factor(name, names(exp_labels), exp_labels))

ggplot(plot_data) +

  
  geom_point(aes(x = median, y = name, colour = factor(.chain)),
             size = 0.8,
             position = position_dodge2(width = 0.3)) +
  geom_linerange(aes(xmin = q95_lower, xmax = q95_upper, y = name, colour = factor(.chain)),
                 position = position_dodge2(width = 0.3)) +
  
  facet_wrap(~variable,
             ncol = 4, scales = "free_x") +
  
  xlab("Value (median, 95% CrI)") + ylab(NULL) +
  
  scale_colour_discrete(type = rep(RColorBrewer::brewer.pal(9, "Blues")[c(6,9)], 10)) +
  
  plot_theme_paper +
  
  theme(strip.text = element_markdown(size = 12, family = "Utopia"),
        panel.grid.major.y = element_gridline,
        axis.text.y = element_markdown(),
        legend.position = "none") +
  
  ggtitle("Ha Nam study inference results")

