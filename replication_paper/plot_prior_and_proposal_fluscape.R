
library(targets)
tar_source()

source("replication_paper/common.R")

summary <- tar_read(combined_summaries) %>%
  filter(exp_group == "prior_proposal_fluscape")

name_order <- c(
  "kucharski_2018", 
  "hay_2024",
  
  "Bernoulli_0.5_uncorrected",
  "Bernoulli_0.5_corrected",
  "BetaBernoulli_1_1_uncorrected",
  "BetaBernoulli_1_1_corrected"
)

name_labels <- c(
  "Kucharski 2018", 
  "Hay 2024",
  "Bernoulli(0.5) (uncorrected)", 
  "Bernoulli(0.5) (corrected)", 
  "BetaBernoulli(1,1) (uncorrected)",
  "BetaBernoulli(1,1) (corrected)"
)


plot_data <- summary %>%
  mutate(name = str_c(prior_description, "_", proposal_name)) %>% 
  filter(name %in% name_order) %>% 
  bind_rows(summaries_previous %>% filter(name == "kucharski_2018")) %>% 
  filter(variable %in% var_names,
         name %in% name_order,
         run_name %in% c("fluscape_2009_neuts", "fluscape_2009_HI")) %>%
  mutate(name = fct_rev(factor(name, name_order, name_labels)),
         variable = factor(variable, var_names, var_labels))

ggplot(plot_data) +
  # annotate("rect", ymin = 0.5, ymax = 2.5,
  #          xmin = -Inf, xmax = Inf, fill = "grey50", alpha = 0.1) +
  # 
  # annotate("rect", ymin = 3.5, ymax = 4.5,
  #          xmin = -Inf, xmax = Inf, fill = ggokabeito::palette_okabe_ito(5), alpha = 0.1) +
  
  geom_point(aes(x = median, y = name),
             size = 0.8,
             position = position_dodge2(width = 0.3)) +
  geom_linerange(aes(xmin = q95_lower, xmax = q95_upper, y = name),
                 position = position_dodge2(width = 0.3)) +
  
  geom_vline(aes(xintercept = median),
             linetype = "14",
             plot_data %>% filter(name == "Kucharski 2018")) +
  
  facet_grid(cols = vars(variable), rows = vars(run_name),
             scales = "free_x") +
  
  xlab("Value (median, 95% CrI)") + ylab(NULL) +
  
  plot_theme_paper +
  
  theme(strip.text = element_markdown(size = 12, family = "Utopia"),
        panel.grid.major.y = element_gridline) +
  
  ggtitle("Fluscape study inference results")
