

library(targets)
tar_source()

source("replication_paper/common.R")

chains <- tar_read(combined_chains)

name_order <- c(
  "Bernoulli_0.5_uncorrected"
  # "Bernoulli_0.5_corrected",
  # "BetaBernoulli_1_1_uncorrected",
  # "BetaBernoulli_1_1_corrected"
)

plot_data <- chains |>
  mutate(name = str_c(prior_description, "_", proposal_name)) |>
  
  filter(exp_group == "prior_proposal_fluscape",
         name %in% name_order) |> 
  
  pivot_longer(any_of(var_names),
               names_to = "variable")

plot_data |>
  ggplot() +
  geom_line(aes(x = .iteration, y = value, colour = factor(.chain))) +
  
  plot_theme_paper +
  
  facet_grid(cols = vars(name), rows = vars(variable), scales = "free_y") +
  
  theme(panel.grid.major = element_gridline)
