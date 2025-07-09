
library(targets)
library(tidyverse)

tar_source()

run_name <- "hanam_2018_age"
chain_name <- "nonlinear_test"

run_dir <- str_c("runs/", run_name, "/")
model_data <- read_model_data(str_c(run_dir, "/model_data.hdf5"))


ppd <- arrow::read_parquet(str_c(run_dir, "ppd_", chain_name, ".parquet")) |>
  process_data_df(model_data$modelled_years)


chain <- read_chain(str_c(run_dir, "chain_", chain_name, ".parquet"))

ix_subject_ex <- 29


subjects <- unique(ppd$ix_subject)

pdf(str_c(run_dir, "ppd_", chain_name, ".pdf"), width = 13, height = 8)
for(ix_subject_ex in subjects) {
  plot_data_ppd <- ppd |> 
    filter(ix_subject == ix_subject_ex,
           ix_sample <= 30)
  
  plot_data_obs <- model_data$observations |>
    filter(ix_subject == ix_subject_ex)
  
  p <- ggplot() +
    
    geom_line(aes(x = strain_year, y = observed_titre, group = ix_sample),
              alpha = 0.5, linewidth = 0.5, colour = "grey30",
              plot_data_ppd) +
    geom_point(aes(x = strain_year, y = observed_titre),
               colour = "red", size = 0.5,
               plot_data_obs) +
    
    scale_x_continuous(labels = function(x) { map_chr(x, function(y) str_c("'", str_sub(y, 3, 4)))}) +
    
    facet_wrap(~year_observed) +
    
    plot_theme_paper +
    theme(axis.line.x = element_line(colour = "grey70"))
  
  plot(p)
}
dev.off()



chain_raw <- arrow::read_parquet(str_c(run_dir, "chain_", chain_name, ".parquet"))

get_lp_mixis(chain_raw, n_warmup = 15000, tar_read(hanam_2018_age), "nonlinear", matrix_beta_bernoulli_1_1, fixed_params = NULL)


