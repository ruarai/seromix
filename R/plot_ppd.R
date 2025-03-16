
source("R/common.R")


complete_obs <- read_parquet("data/complete_obs.parquet") %>%
  filter(i < 5)

obs_df <- read_parquet("data/obs.parquet") %>%
  filter(i < 5)

ggplot() +
  geom_line(aes(x = t, y = y, group = s, colour = factor(s)),
            complete_obs) +
  geom_point(aes(x = t, y = y, group = s, colour = factor(s)),
             obs_df) +
  facet_wrap(~i, ncol = 1) +
  
  scale_x_continuous(breaks = 1:10)



chain_df <- read_draws("data/chain.parquet") %>%
  filter()

warmup_steps <- 43000
thinning <- 10

mcmc_trace(chain_df, "mu_sum", n_warmup = warmup_steps / thinning)

parnames <- colnames(chain_df)[5:10]

chain_df_filt <- chain_df %>%
  filter(.iteration > warmup_steps)


mcmc_hist(chain_df_filt, parnames)
mcmc_pairs(chain_df_filt,
           parnames[1:5],
           off_diag_args = list(size = 0.5))


infections_df <- read_parquet("data/infections_df.parquet") %>%
  `colnames<-`(c("t", "i"))

plot_data_inf <- chain_df_filt %>%
  spread_draws(infections[t, i]) %>%
  group_by(t, i, .chain) %>%
  summarise(p = 1 - sum(infections == 0) / n()) %>%
  left_join(infections_df %>% mutate(inf = TRUE)) %>%
  mutate(inf = replace_na(inf, FALSE))

autoplot(precrec::evalmod(scores = plot_data_inf$p, labels = plot_data_inf$inf))

plot_data_inf %>%
  filter(i < 50) %>% 
  ggplot() +
  geom_tile(aes(x = t, y = i, fill = p)) +
  
  geom_point(aes(x = t, y = i), infections_df %>% filter(i < 50), colour = "red")





ppd_df <- read_draws("data/ppd.parquet")
ppd_obs <- read_parquet("data/ppd_obs.parquet")

plot_data <- ppd_df %>%
  filter(.draw %% 10 == 0) %>%
  spread_draws(obs_titre[ix_row]) %>%
  left_join(ppd_obs %>% mutate(ix_row = row_number())) %>%
  filter(i < 3)


ggplot() +
  geom_line(aes(x = t, y = obs_titre, group = .draw),
            alpha = 0.3,
            plot_data) +
  
  geom_point(aes(x = t, y = y),
             obs_df %>% filter(i < 3), colour = "red") +
  
  geom_vline(aes(xintercept = t), infections_df %>% filter(i < 2) %>% mutate(s = t),
             colour = "red") +
  
  facet_wrap(~i * s, ncol = 8)

