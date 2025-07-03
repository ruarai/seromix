
library(targets)
tar_source()

source("replication_paper/common.R")

summary <- tar_read(combined_summaries) |> 
  filter(exp_group == "titre_correction")

plot_data <- summary |>
  mutate(name = str_c(prior_description, "_", proposal_name)) |> 
  filter(initial_params_name == "kucharski_data_study",
         prior_description == "BetaBernoulli_1_1",
         proposal_name == "corrected") |> 
  
  bind_rows(summaries_previous |> filter(name == "kucharski_2018")) |> 
  filter(run_name == "hanam_2018",
         variable %in% var_names,
         name %in% name_order) |>
  mutate(name = factor(if_else(use_corrected_titre, "Corrected", "Uncorrected")),
         variable = factor(variable, var_names, var_labels))

ggplot(plot_data) +
  annotate("rect", ymin = 0.5, ymax = 1.5,
           xmin = -Inf, xmax = Inf, fill = "grey50", alpha = 0.1) +
  
  geom_point(aes(x = median, y = name),
             size = 0.8,
             position = position_dodge2(width = 0.3)) +
  geom_linerange(aes(xmin = q95_lower, xmax = q95_upper, y = name),
                 position = position_dodge2(width = 0.3)) +
  
  facet_wrap(~variable,
             ncol = 4, scales = "free_x") +
  
  xlab("Value (median, 95% CrI)") + ylab(NULL) +
  
  # coord_cartesian(xlim = c(0, NA)) +
  
  plot_theme_paper +
  
  theme(strip.text = element_markdown(size = 12, family = "Utopia"),
        panel.grid.major.y = element_gridline)


ggsave(
  "replication_paper/results/corrected_titre.png",
  device = png,
  width = 12,
  height = 4, bg = "white"
)



