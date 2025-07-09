
library(targets)
tar_source()

source("replication_paper/common.R")

summary <- tar_read(combined_summaries) |>
  filter(exp_group == "model_comparison", mixture_importance_sampling)

lp_mixis_summary <- tar_read(combined_lp_mixis) %>%
  drop_na(lp_mixis) %>%
  mutate(variable = "mixis_lp", median = lp_mixis)

exp_labels <- c(
  "with_both" = "With both",
  "with_intercept" = "With intercept",
  "with_age_effect" = "With age-effect",
  "kucharski" = "Baseline"#,
  # "without_seniority" = "Without seniority"
)

plot_data <- summary |>
  
  bind_rows(lp_mixis_summary) |> 
  mutate(name = exp_name) |>
  filter(name %in% names(exp_labels)) |> 
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
        panel.background = element_facet_background,
        legend.position = "none") +
  
  ggtitle("Ha Nam study inference results")

