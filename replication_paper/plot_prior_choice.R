
library(targets)
tar_source()

source("replication_paper/common.R")

summary <- tar_read(combined_summaries) |>
  filter(exp_group == "prior_comparison")

name_order <- c(
  "BetaBernoulli_1_1",
  "BetaBernoulli_2.5_8",
  "BetaBernoulliSubjectVarying_1_1",
  "BetaBernoulliSubjectVarying_2.5_8",
  "BetaBernoulliTimeVarying_1_1",
  "BetaBernoulliTimeVarying_2.5_8",
  "Bernoulli_0.5",
  "Bernoulli_0.3",
  "Bernoulli_0.1"
)

name_labels <- c(
  "BetaBernoulli(1, 1)",
  "BetaBernoulli(2.5, 8)",
  "BetaBernoulli<sub>SV</sub>(1, 1)",
  "BetaBernoulli<sub>SV</sub>(2.5, 8)",
  "BetaBernoulli<sub>TV</sub>(1, 1)",
  "BetaBernoulli<sub>TV</sub>(2.5, 8)",
  "Bernoulli(0.5)", 
  "Bernoulli(0.3)", 
  "Bernoulli(0.1)"
)

plot_data <- summary |>
  mutate(name = prior_description) |> 
  filter(run_name == "hanam_2018",
         variable %in% names(var_labels)) |>
  mutate(name = fct_rev(factor(name, name_order, name_labels)),
         variable = factor(variable, names(var_labels), var_labels))

ggplot(plot_data) +
  
  geom_rect(aes(xmin = -Inf, xmax = Inf, ymin = y - 0.5, ymax = y + 0.5),
            tibble(y = seq(1, length(name_labels), by = 2)), fill = "grey50", alpha = 0.1) +
  
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



ggsave(
  "replication_paper/results/prior_choice.png",
  device = png,
  width = 12,
  height = 8, bg = "white"
)


