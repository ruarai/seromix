
source("R/common.R")




chain_df <- read_draws("data/chain.parquet")



ppd_df <- read_draws("data/ppd_df.parquet")



plot_data <- ppd_df %>%
  spread_draws(titer_obs[t])


ggplot() +
  geom_point(aes(x = t, y = titer_obs, group = .draw),
             plot_data) 
