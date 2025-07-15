
library(targets)
tar_source()

source("replication_paper/common.R")


tar_read(combined_singular_summaries) |>
  filter(exp_group == "prior_proposal") |>
  filter(rhat > 1.1) |> 
  select(variable, prior_description, proposal_name, rhat)

summary <- tar_read(combined_summaries) |>
  filter(exp_group == "prior_proposal")

name_labels <- c(
  "kucharski_2018" = "Kucharski 2018",
  # "hay_2024" = "Hay 2024",
  
  "Bernoulli_0.5_uncorrected" = "Bernoulli(0.5), uncorrected",
  "Bernoulli_0.5_corrected" = "Bernoulli(0.5), corrected",
  "BetaBernoulli_1_1_uncorrected" = "BetaBernoulli(1,1), uncorrected",
  "BetaBernoulli_1_1_corrected" = "BetaBernoulli(1,1), corrected",
  "pigeons" = "BetaBernoulli(1,1), corrected (PT)"
)

summary_pigeons <- bind_rows(
  arrow::read_parquet("runs/hanam_2018/chain_pigeons_6_1.parquet") |> mutate(chain = 1),
  arrow::read_parquet("runs/hanam_2018/chain_pigeons_6_2.parquet") |> mutate(chain = 2),
) |> 
  reformat_pigeons_chain(tar_read(hanam_2018)) |> 
  summarise_chain(0, tar_read(hanam_2018), by_chain = TRUE) |> 
  mutate(run_name = "hanam_2018",
         name = "pigeons")


plot_data <- summary |>
  mutate(name = str_c(prior_description, "_", proposal_name)) |> 
  bind_rows(summaries_previous |> filter(name == "kucharski_2018") |> mutate(.chain = 1),
            summary_pigeons) |> 
  filter(run_name == "hanam_2018",
         variable %in% names(var_labels),
         name %in% names(name_labels)) |>
  mutate(name = fct_rev(factor(name, names(name_labels), name_labels)),
         variable = factor(variable, names(var_labels), var_labels))




ggplot(plot_data) +
  annotate("rect", ymin = 0.5, ymax = 1.5,
           xmin = -Inf, xmax = Inf, fill = ggokabeito::palette_okabe_ito(6), alpha = 0.1) +
  annotate("rect", ymin = 1.5, ymax = 2.5,
           xmin = -Inf, xmax = Inf, fill = "grey50", alpha = 0.1) +
  annotate("rect", ymin = 3.5, ymax = 4.5,
           xmin = -Inf, xmax = Inf, fill = "grey50", alpha = 0.1) +
  annotate("rect", ymin = 5.5, ymax = 6.5,
           xmin = -Inf, xmax = Inf, fill = ggokabeito::palette_okabe_ito(5), alpha = 0.1) +
  geom_linerange(aes(xmin = q95_lower, xmax = q95_upper, y = name, colour = factor(.chain)),
                 position = position_dodge2(width = 0.3)) +
  
  geom_point(aes(x = median, y = name, colour = factor(.chain)),
             size = 0.8,
             position = position_dodge2(width = 0.3)) +
  
  geom_vline(aes(xintercept = median),
             linetype = "14",
             plot_data |> filter(name == "Kucharski 2018")) +
  
  facet_wrap(~variable,
             ncol = 4, scales = "free_x") +
  
  xlab("Value (median, 95% CrI)") + ylab(NULL) +

  scale_colour_discrete(type = rep(RColorBrewer::brewer.pal(9, "Blues")[c(6,9)], 10)) +
  
  
  plot_theme_paper +
  
  theme(strip.text = element_markdown(size = 12, family = "Utopia"),
        panel.grid.major.y = element_gridline,
        legend.position = "none") +
  
  ggtitle("Ha Nam study inference results")


ggsave(
  "replication_paper/results/prior_and_proposal.png",
  device = png,
  width = 12,
  height = 8, bg = "white"
)



