

tar_source() 
run_name <- "fluscape_2009_neuts"
run_dir <- str_c("runs/", run_name)

chain_files <- list.files(path = run_dir, pattern = "chain_prior_*", full.names = TRUE)


col_vars <- c("mu_long", "mu_short", "omega", "sigma_long", "sigma_short", "tau", "obs_sd")

chain_data <- chain_files |>
  map(
    function(x) {
      x |> 
        read_chain() |>
        filter(.iteration > 40000) |>
        select(any_of(c(".iteration", ".chain", col_vars)))
    }
  )

chain_name_order <- c("prior_50_uncorrected", "prior_15_uncorrected",
                      "prior_50_corrected", "prior_15_corrected",
                      "prior_beta_1.3_8.0_corrected")

reported_values <- tribble(
  ~name, ~value, ~run_name,
  "mu_long", 2.02, "hanam_2018",
  "mu_short", 2.69, "hanam_2018",
  "sigma_long", 0.13, "hanam_2018",
  "sigma_short", 0.03, "hanam_2018",
  "tau", 0.039, "hanam_2018",
  "omega", 0.79, "hanam_2018",
  "obs_sd", 1.29, "hanam_2018",
  
  "mu_long", 0.97, "fluscape_2009_HI",
  "sigma_long", 0.099, "fluscape_2009_HI",
  "tau", 0.016, "fluscape_2009_HI",
  "obs_sd", 1.5, "fluscape_2009_HI",
  
  "mu_long", 1.38, "fluscape_2009_neuts",
  "sigma_long", 0.130, "fluscape_2009_neuts",
  "tau", 0.02, "fluscape_2009_neuts",
  "obs_sd", 1.69, "fluscape_2009_neuts"
) |>
  filter(run_name == !!run_name)


plot_data_chain <- chain_data |>
  bind_rows(.id = "chain_file") |>
  mutate(chain_file = chain_files[as.numeric(chain_file)],
         chain_name = str_match(chain_file, "chain_(.+)\\.parquet")[,2],
         chain_name = fct_rev(factor(chain_name, chain_name_order))) |>
  pivot_longer(any_of(col_vars))

if(str_starts(run_name, "fluscape_2009")) {
  plot_data_chain <- plot_data_chain |>
    filter(name != "sigma_short")
}


ggplot() +
  ggdist::stat_interval(aes(x = value, y = chain_name),
                        position = position_dodge(), size = 3.0,
                        .width = c(0.5, 0.75, 0.9, 0.95, 0.99),
                        plot_data_chain) +
  
  geom_vline(aes(xintercept = value), linetype = "44",
             reported_values) +
  
  scale_colour_brewer() +
  
  facet_wrap(~name, scales = "free_x")

