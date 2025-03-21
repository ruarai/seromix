
source("R/common.R")

true_values <- tribble(
  ~name, ~value,
  "mu_long", 2.0,
  "mu_sum", 2.0 + 2.0,
  "mu_short", 2.0,
  "sigma_long", 0.2,
  "sigma_short", 0.1,
  "tau", 0.05
)

run_dir <- "runs/sim_study_hanam_2018_2/"

infections_df <- read_parquet(str_c(run_dir, "infections_df.parquet")) %>%
  rename(ix_t = t)

subject_birth_data <- rhdf5::h5read(str_c(run_dir, "model_data.hdf5"), "subject_birth_data")
modelled_years <- rhdf5::h5read(str_c(run_dir, "model_data.hdf5"), "modelled_years")

observations_df <- rhdf5::h5read(str_c(run_dir, "model_data.hdf5"), "observations_df")
complete_obs <- rhdf5::h5read(str_c(run_dir, "model_data.hdf5"), "complete_obs")

chain_df <- read_draws(str_c(run_dir, "chain.parquet"))

warmup_steps <- 18000

parnames <- c(colnames(chain_df)[5:9], "mu_short")

chain_df %>%
  select(.iteration, .chain, any_of(parnames)) %>%
  pivot_longer(any_of(parnames)) %>%
  ggplot() +
  geom_line(aes(x = .iteration, y = value, colour = factor(.chain))) +
  
  geom_vline(xintercept = warmup_steps) +
  
  facet_wrap(~name, scales = "free_y")



chain_inf_draws <- chain_df %>%
  
  spread_draws(infections[ix_t, ix_subject]) %>%
  
  left_join(subject_birth_data) %>%
  filter(ix_t > ix_t_birth)


plot_data_inf_chain <- chain_inf_draws %>%
  group_by(.iteration, .chain) %>%
  
  summarise(n_inf = sum(infections))

ggplot() +
  geom_line(aes(x = .iteration, y = n_inf, group = .chain, colour = factor(.chain)),
            plot_data_inf_chain) +
  
  geom_hline(yintercept = nrow(infections_df))


chain_df_filt <- chain_df %>%
  filter(.iteration > warmup_steps) %>%
  mutate(mu_short = mu_sum - mu_long)

chain_filt_wide <- chain_df_filt %>% 
  select(.iteration, .chain, any_of(parnames)) %>%
  pivot_longer(any_of(parnames))

mcmc_pairs(chain_df_filt, parnames)

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
  
  group_by(ix_t, ix_subject) %>%
  summarise(p = 1 - sum(infections == 0) / n()) %>%
  left_join(infections_df %>% mutate(inf = TRUE)) %>%
  mutate(inf = replace_na(inf, FALSE)) %>%
  
  left_join(subject_birth_data) %>%
  filter(ix_t > ix_t_birth)



plot_data_inf %>%
  # filter(ix_subject < 50) %>% 
  ggplot() +
  geom_tile(aes(x = ix_t, y = ix_subject, fill = p)) +
  
  geom_point(aes(x = ix_t, y = ix_subject), infections_df, colour = "red")


autoplot(precrec::evalmod(scores = plot_data_inf$p, labels = plot_data_inf$inf), curvetype = "ROC")



plot_data_prevalence <- chain_df_filt %>%
  
  spread_draws(infections[ix_t, ix_subject]) %>%
  
  group_by(ix_t) %>%
  summarise(p = 1 - sum(infections == 0) / n()) %>%
  
  mutate(year = modelled_years[ix_t]) %>% 
  
  left_join(subjects_alive)


ggplot() +
  geom_line(aes(x = year, y = p / p_alive),
            plot_data_prevalence)


ix_subject_ex <- 45

ppd_obs <- read_parquet(str_c(run_dir, "ppd_obs.parquet")) %>%
  filter(ix_subject == ix_subject_ex) %>%
  mutate(year_sampled = modelled_years[ix_t_obs],
         strain_year = modelled_years[ix_strain]) %>%
  
  filter(year_sampled %% 2 == 0)


ggplot() +
  
  geom_point(aes(x = strain_year, y = observed_titre),
             size = 0.5,
             observations_df %>% filter(ix_subject == ix_subject_ex)) +

  geom_line(aes(x = strain_year, y = observed_titre, group = draw),
            linewidth = 0.3, alpha = 0.3,
            ppd_obs) +
  
  # geom_line(aes(x = strain_year, y = observed_titre),
  #           colour = "red", alpha = 1.0,
  #           complete_obs %>% filter(ix_subject == ix_subject_ex)) +
  
  
  facet_wrap(~year_sampled, ncol = 3)
