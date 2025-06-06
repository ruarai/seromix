
library(targets)
library(tidyverse)

tar_source()

summary <- tar_read(combined_summaries)

var_names <- c("mu_long", "mu_short", "sigma_long", "sigma_short",  "tau", "omega", "obs_sd", "total_inf")
var_labels <- c(
  "<i>μ</i><sub>long</sub> (long-term boost)", 
  "<i>μ</i><sub>short</sub> (short-term boost)",
  "<i>σ</i><sub>long</sub> (short-term cross.)", 
  "<i>σ</i><sub>short</sub> (long-term cross.)",
  "<i>τ</i> (seniority)",
  "<i>ω</i> (waning)",
  "obs sd.", 
  "Total inf."
)

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
  bind_rows(summaries_previous %>% filter(name == "kucharski_2018")) %>% 
  filter(run_name == "hanam_2018",
         variable %in% var_names) %>%
  mutate(name = fct_rev(factor(name, name_order, name_labels)),
         variable = factor(variable, var_names, var_labels))



ggplot(plot_data) +
  annotate("rect", ymin = 0.5, ymax = 1.5,
           xmin = -Inf, xmax = Inf, fill = "grey50", alpha = 0.1) +
  annotate("rect", ymin = 2.5, ymax = 3.5,
           xmin = -Inf, xmax = Inf, fill = "grey50", alpha = 0.1) +
  annotate("rect", ymin = 4.5, ymax = 5.5,
           xmin = -Inf, xmax = Inf, fill = ggokabeito::palette_okabe_ito(5), alpha = 0.1) +
  
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
        panel.grid.major.y = element_gridline)


df <- tar_read(chain_hanam_2018_uncorrected_Bernoulli_0.5) %>%
  clean_chain()

ggplot(df) +
  geom_line(aes(x = .iteration, y = mu_long, colour = factor(.chain)))

