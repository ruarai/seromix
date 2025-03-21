
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

run_dir <- "runs/sim_study_simple_1/"

model_data <- tar_read(sim_study_1_data)
fit_data <- read_fit_data(str_c(run_dir, "fit_data.hdf5"), model_data$modelled_years)

chain <- fit_data$chain
ppd <- fit_data$ppd

chain <- chain %>%
  mutate(mu_short = mu_sum - mu_long)

warmup_steps <- 1800

parnames <- c(colnames(chain)[5:9], "mu_short")

chain %>%
  select(.iteration, .chain, any_of(parnames)) %>%
  pivot_longer(any_of(parnames)) %>%
  ggplot() +
  geom_line(aes(x = .iteration, y = value, colour = factor(.chain))) +
  
  geom_vline(xintercept = warmup_steps) +
  
  facet_wrap(~name, scales = "free_y")

chain_filt <- chain %>%
  filter(.iteration > warmup_steps)

chain_filt_wide <- chain_filt %>% 
  select(.iteration, .chain, any_of(parnames)) %>%
  pivot_longer(any_of(parnames))

ggplot() +
  geom_line(aes(x = .iteration, y = value, colour = factor(.chain)),
            chain_filt_wide) +
  facet_wrap(~name, scales = "free_y")

ggplot() +
  geom_histogram(aes(x = value), bins = 50,
                 chain_filt_wide) +
  
  geom_vline(aes(xintercept = value), true_values) +
  
  facet_wrap(~name, scales = "free")

plot_data_inf <- chain_filt %>%
  
  spread_draws(infections[ix_t, ix_subject]) %>%
  
  group_by(ix_t, ix_subject) %>%
  summarise(p = 1 - sum(infections == 0) / n()) %>%
  left_join(model_data$infections %>% mutate(inf = TRUE)) %>%
  mutate(inf = replace_na(inf, FALSE)) %>%
  
  left_join(model_data$subject_birth_data) %>%
  mutate(ix_t_birth = replace_na(ix_t_birth, 0)) %>% 
  filter(ix_t >= ix_t_birth)



plot_data_inf %>%
  ggplot() +
  geom_tile(aes(x = ix_t, y = ix_subject, fill = p)) +
  
  geom_point(aes(x = ix_t, y = ix_subject), model_data$infections, colour = "red")


autoplot(precrec::evalmod(scores = plot_data_inf$p, labels = plot_data_inf$inf), curvetype = "ROC")


ix_subject_ex <- 15

