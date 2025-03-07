
source("R/common.R")



obs_df <- read_parquet("data/obs.parquet")

chain_df <- read_draws("data/chain.parquet")

parnames <- colnames(chain_df)[4:8]

bayesplot::mcmc_hist(chain_df, parnames)



ppd_df <- read_draws("data/ppd.parquet")
ppd_obs_df <- read_parquet("data/ppd_obs.parquet")

plot_data <- ppd_df %>%
  filter(.draw < 20) %>%
  spread_draws(obs_titre[ix_row]) %>%
  left_join(ppd_obs_df %>% mutate(ix_row = row_number()))

ggplot() +
  geom_tile(aes(x = t, y = ix_strain, fill = obs_titre),
            plot_data %>% filter(.draw == 1))


ggplot() +
  geom_line(aes(x = t, y = obs_titre, group = .draw),
            plot_data) +
  
  geom_point(aes(x = t, y = titre),
             obs_df, colour = "red") +
  
  facet_wrap(~ix_strain)
