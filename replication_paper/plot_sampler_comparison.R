
library(targets)
tar_source()

source("replication_paper/common.R")

summary <- tar_read(combined_summaries) |>
  filter(exp_group == "sampler_comparison")


plot_data <- summary |>
  
  mutate(name = sampler_name) |>
  mutate(variable = factor(variable, names(var_labels), var_labels))

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
