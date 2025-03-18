
source("R/common.R")

true_values <- tribble(
  ~name, ~value,
  "mu_long", 2.0,
  "mu_sum", 2.0 + 2.7,
  "sigma_long", 0.13,
  "sigma_short", 0.03,
  "tau", 0.04
)

run_dir <- "runs/sim_inf_hanam_2018/"


complete_obs <- read_parquet(str_c(run_dir, "complete_obs.parquet")) %>%
  mutate(year_sampled = modelled_years[ix_t_obs],
         strain_year = modelled_years[ix_strain]) %>%
  
  filter(year_sampled %% 2 == 0)

obs_df <- read_parquet(str_c(run_dir, "obs.parquet")) %>%
  mutate(year_sampled = modelled_years[ix_t_obs],
         strain_year = modelled_years[ix_strain])

infections_df <- read_parquet(str_c(run_dir, "infections_df.parquet")) %>% 
  `colnames<-`(c("ix_t", "ix_subject")) %>%
  mutate(strain_year = modelled_years[ix_t])

ggplot() +
  geom_line(aes(x = strain_year, y = observed_titre),
                 complete_obs %>% filter(ix_subject == 50)) +
  
  geom_point(aes(x = strain_year, y = observed_titre),
             size = 0.5,
             obs_df %>% filter(ix_subject == 50)) +
  
  xlab("Strain year") +
  
  facet_wrap(~year_sampled, ncol = 3)


chain_df <- read_draws(str_c(run_dir, "chain.parquet"))

warmup_steps <- 3000
thinning <- 1

parnames <- colnames(chain_df)[5:9]

chain_df %>%
  select(.iteration, .chain, all_of(parnames)) %>%
  pivot_longer(all_of(parnames)) %>%
  ggplot() +
  geom_line(aes(x = .iteration, y = value, colour = factor(.chain))) +
  
  geom_vline(xintercept = warmup_steps) +
  
  facet_wrap(~name, scales = "free_y")

chain_df_filt <- chain_df %>%
  filter(.iteration > warmup_steps)

chain_filt_wide <- chain_df_filt %>%
  select(.iteration, .chain, all_of(parnames)) %>%
  pivot_longer(all_of(parnames))

ggplot() +
  geom_line(aes(x = .iteration, y = value, colour = factor(.chain)),
            chain_filt_wide) +
  facet_wrap(~name, scales = "free_y")

ggplot() +
  geom_histogram(aes(x = value), bins = 50,
                 chain_filt_wide) +
  
  geom_vline(aes(xintercept = value), true_values) +
  
  facet_wrap(~name, scales = "free")



plot_data_inf <- chain_df_filt %>%
  
  spread_draws(infections[ix_t, ix_subject]) %>%
  group_by(ix_t, ix_subject, .chain) %>%
  summarise(p = 1 - sum(infections == 0) / n()) %>%
  left_join(infections_df %>% mutate(inf = TRUE)) %>%
  mutate(inf = replace_na(inf, FALSE)) %>%
  
  left_join(subject_birth_data) %>%
  filter(ix_t > ix_t_birth)

autoplot(precrec::evalmod(scores = plot_data_inf$p, labels = plot_data_inf$inf))

plot_data_inf %>%
  filter(ix_subject < 50) %>% 
  ggplot() +
  geom_tile(aes(x = ix_t, y = ix_subject, fill = p)) +
  
  geom_point(aes(x = ix_t, y = ix_subject), infections_df %>% filter(ix_subject < 50), colour = "red")




ix_subject_ex <- 45

ppd_obs <- read_parquet(str_c(run_dir, "ppd_obs.parquet")) %>%
  filter(ix_subject == ix_subject_ex) %>%
  mutate(year_sampled = modelled_years[ix_t_obs],
         strain_year = modelled_years[ix_strain]) %>%
  
  filter(year_sampled %% 2 == 0)


ggplot() +
  
  geom_point(aes(x = strain_year, y = observed_titre),
             size = 0.5,
             obs_df %>% filter(ix_subject == ix_subject_ex)) +

  geom_line(aes(x = strain_year, y = observed_titre, group = draw),
            linewidth = 0.3, alpha = 0.3,
            ppd_obs) +
  
  geom_line(aes(x = strain_year, y = observed_titre),
            colour = "red", alpha = 1.0,
            complete_obs %>% filter(ix_subject == ix_subject_ex)) +
  
  
  facet_wrap(~year_sampled, ncol = 3)
