
source("R/common.R")


complete_obs <- read_parquet("data/complete_obs.parquet") %>%
  filter(i < 10)

obs_df <- read_parquet("data/obs.parquet") %>%
  filter(i < 10)

ggplot() +
  geom_line(aes(x = t, y = y, group = s, colour = factor(s)),
            complete_obs) +
  geom_point(aes(x = t, y = y, group = s, colour = factor(s)),
             obs_df) +
  facet_wrap(~i, ncol = 1) +
  
  scale_x_continuous(breaks = 1:10)



chain_df <- read_draws("data/chain.parquet")

parnames <- colnames(chain_df)[4:9]


bayesplot::mcmc_hist(chain_df, parnames)
bayesplot::mcmc_pairs(chain_df %>% sample_n(500),
                      parnames[1:5],
                      off_diag_args = list(size = 0.5))


infections_df <- read_parquet("data/infections_df.parquet") %>%
  `colnames<-`(c("t", "i"))

chain_df %>%
  spread_draws(infections[t, i]) %>%
  filter(i < 10) %>%
  group_by(t, i, .chain) %>%
  summarise(p = 1 - sum(infections == 0) / n()) %>%
  ggplot() +
  geom_line(aes(x = t, y = p, group = .chain, colour = factor(.chain))) +
  
  geom_vline(aes(xintercept = t), infections_df %>% filter(i < 10)) +
  facet_wrap(~i, ncol = 1) +
  scale_x_continuous(breaks = 1:10)


chain_df %>%
  spread_draws(infections[t, i]) %>%
  group_by(t, i) %>%
  summarise(p = 1 - sum(infections == 0) / n()) %>%
  ggplot() +
  geom_tile(aes(x = t, y = i, fill = p)) +
  
  geom_point(aes(x = t, y = i), infections_df, colour = "red")


ppd_df <- read_draws("data/ppd.parquet")
ppd_obs <- read_parquet("data/ppd_obs.parquet")

plot_data <- ppd_df %>%
  filter(.draw %% 10 == 0) %>%
  spread_draws(obs_titre[ix_row]) %>%
  left_join(ppd_obs %>% mutate(ix_row = row_number())) %>%
  filter(i < 4)


ggplot() +
  geom_line(aes(x = t, y = obs_titre, group = .draw),
            alpha = 0.3,
            plot_data) +
  
  geom_point(aes(x = t, y = y),
             obs_df %>% filter(i < 4), colour = "red") +
  
  geom_vline(aes(xintercept = t), infections_df %>% filter(i < 4) %>% mutate(s = t),
             colour = "red") +
  
  facet_wrap(~i * s, ncol = 6)

