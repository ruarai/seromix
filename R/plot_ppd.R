
source("R/common.R")



obs_df <- read_parquet("data/obs.parquet")

chain_df <- read_draws("data/chain.parquet")

parnames <- colnames(chain_df)[4:9]


bayesplot::mcmc_hist(chain_df, parnames)
# bayesplot::mcmc_pairs(chain_df, parnames)


chain_df %>%
  filter(.iteration > 4000) %>% 
  spread_draws(infections[t, ix_strain]) %>%
  group_by(t, ix_strain) %>%
  summarise(p = 1 - sum(infections == 0) / n()) %>%
  ggplot() +
  geom_line(aes(x = t, y = p, group = ix_strain))


ppd_df <- read_draws("data/ppd.parquet")
ppd_obs_df <- read_parquet("data/ppd_obs.parquet")

plot_data <- ppd_df %>%
  filter(.draw %% 100 == 0) %>%
  spread_draws(obs_titre[ix_row]) %>%
  left_join(ppd_obs_df %>% mutate(ix_row = row_number()))


ggplot() +
  geom_line(aes(x = t, y = obs_titre, group = .draw),
            alpha = 0.1,
            plot_data) +
  
  geom_point(aes(x = t, y = titre),
             obs_df, colour = "red") +
  
  facet_wrap(~ix_strain)

