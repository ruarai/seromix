
library(targets)
tar_source()

source("replication_paper/common.R")

summary <- tar_read(combined_summaries) |>
  filter(exp_group == "model_comparison")

plot_data <- summary |>
  mutate(name = initial_params_name)

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
