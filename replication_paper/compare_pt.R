
library(targets)
tar_source()

source("replication_paper/common.R")



summary_pigeons <- arrow::read_parquet("runs/hanam_2018/chain_pigeons_5_mixis.parquet") |> 
  reformat_pigeons_chain(tar_read(hanam_2018)) |> 
  summarise_chain(0, tar_read(hanam_2018), by_chain = FALSE) |> 
  bind_rows(tibble(variable = "mixis_lp", median = -16791.027416869154)) |> 
  mutate(run_name = "hanam_2018_pt",
         name = "pigeons",
         .chain = 1,
         variable = if_else(variable == "log_density", "lp", variable))

summary <- tar_read(combined_summaries) |>
  filter(exp_group == "baseline_mixis", run_name != "hanam_2018_age")

lp_mixis_summary <- tar_read(combined_lp_mixis) %>%
  filter(exp_group == "baseline_mixis", run_name != "hanam_2018_age") |> 
  mutate(variable = "mixis_lp", median = lp_mixis)


plot_data <- summary |>
  
  bind_rows(lp_mixis_summary) |> 
  bind_rows(summary_pigeons) |> 
  mutate(variable = factor(variable, names(var_labels), var_labels))


ggplot(plot_data) +
  geom_point(aes(x = median, y = run_name, colour = factor(.chain)),
             size = 0.8,
             position = position_dodge2(width = 0.3)) +
  
  geom_linerange(aes(xmin = q95_lower, xmax = q95_upper, y = run_name, colour = factor(.chain)),
                 position = position_dodge2(width = 0.3)) +
  
  
  facet_wrap(~variable,
             ncol = 4, scales = "free_x") +
  
  scale_colour_discrete(type = rep(RColorBrewer::brewer.pal(9, "Blues")[c(6,9)], 10)) +
  
  plot_theme_paper +
  
  theme(strip.text = element_markdown(size = 12, family = "Utopia"),
        panel.grid.major.y = element_gridline,
        axis.text.y = element_markdown(),
        panel.background = element_facet_background,
        legend.position = "none")
