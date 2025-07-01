
library(targets)
tar_source()

source("replication_paper/common.R")



chain_pigeon <- arrow::read_parquet("runs/hanam_2018/chain_pigeons_3.parquet")

model_data <- tar_read(hanam_2018)

n_subjects <- length(model_data$age_distribution)
n_t_steps <- length(model_data$modelled_years)

inf_columns_new <- colnames(chain_pigeon) |> 
  keep(~ str_detect(.x, "infections")) |> 
  map_dbl(~ as.numeric(str_extract(.x, "\\d+")) - 1) |> 
  map_chr(~ str_c("infections[", .x %% n_t_steps + 1, ",", .x %/% n_t_steps + 1, "]"))

col_names <- colnames(chain_pigeon)
col_names[str_detect(col_names, "infections")] <- inf_columns_new

colnames(chain_pigeon) <- col_names


summ_pigeon <- chain_pigeon |> 
  # select(-starts_with("infections")) |> 
  summarise_chain(0, tar_read(hanam_2018), by_chain = FALSE) |> 
  mutate(run_name = "hanam_2018",
         name = "pigeons")


summary <- tar_read(combined_summaries) |>
  filter(exp_group == "prior_proposal")




plot_data <- summary |>
  # mutate(name = str_c(prior_description, "_", proposal_name)) |> 
  bind_rows(summaries_previous |> filter(name == "kucharski_2018")) |> 
  bind_rows(summ_pigeon) |> 
  filter(run_name == "hanam_2018",
         variable %in% var_names,
         # name %in% name_order
         ) |>
  mutate(#name = fct_rev(factor(name, name_order, name_labels)),
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
  
  geom_vline(aes(xintercept = median),
             linetype = "14",
             plot_data |> filter(name == "Kucharski 2018")) +
  
  facet_wrap(~variable,
             ncol = 4, scales = "free_x") +
  
  xlab("Value (median, 95% CrI)") + ylab(NULL) +
  
  plot_theme_paper +
  
  theme(strip.text = element_markdown(size = 12, family = "Utopia"),
        panel.grid.major.y = element_gridline) +
  
  ggtitle("Ha Nam study inference results")
