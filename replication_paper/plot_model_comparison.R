
library(targets)
tar_source()

source("replication_paper/common.R")

summary <- tar_read(combined_summaries) |>
  filter(exp_group == "model_comparison", !mixture_importance_sampling)

lp_mixis_summary <- tar_read(combined_lp_mixis) %>%
  drop_na(lp_mixis) %>%
  mutate(variable = "lp_mixis", median = lp_mixis)

model_names <- c(
  "model_comparison_hanam_2018_age_1" = "no_tau",
  "model_comparison_hanam_2018_age_2" = "no_tau", 
  "model_comparison_hanam_2018_age_3" = "kucharski", 
  "model_comparison_hanam_2018_age_4" = "kucharski",
  "model_comparison_hanam_2018_age_5" = "age_effect",
  "model_comparison_hanam_2018_age_6" = "age_effect"
)

plot_data <- summary |>
  bind_rows(lp_mixis_summary) |> 
  mutate(name = model_names[name])

ggplot(plot_data) +

  
  geom_point(aes(x = median, y = name),
             size = 0.8,
             position = position_dodge2(width = 0.3)) +
  geom_linerange(aes(xmin = q95_lower, xmax = q95_upper, y = name),
                 position = position_dodge2(width = 0.3)) +
  
  facet_wrap(~variable,
             ncol = 4, scales = "free_x") +
  
  xlab("Value (median, 95% CrI)") + ylab(NULL) +
  
  plot_theme_paper +
  
  theme(strip.text = element_markdown(size = 12, family = "Utopia"),
        panel.grid.major.y = element_gridline,
        axis.text.y = element_markdown(),
        legend.position = "none") +
  
  ggtitle("Ha Nam study inference results")
