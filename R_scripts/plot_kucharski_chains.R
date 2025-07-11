
run_data <- read_delim("test2.txt", ";", col_names = c(".iteration", "pr_inf", "pr_param", "eps", "t", "t2")) |> 
  mutate(pr_param = as.numeric(pr_param),
         eps = as.numeric(eps),
         .iteration = as.numeric(.iteration) * 2) 

e <- new.env()

run_data_orig <- map(
  1:8,
  function(ix_chain) {
    load(
      str_c("../flu-model-kucharski/sero_model/posterior_sero_runs/outputR_f2007_2008_2009_2010_2011_2012_s", 
            ix_chain, "_H3HN_linTRUE.RData")
      ,
      envir = e
    )
    tibble(.iteration = (1:length(e$accept_rate_tab)), pr_param = e$accept_rate_tab, eps = e$eps_tab, .chain = ix_chain)
  }
) |> bind_rows()


ggplot() +
  
  geom_step(aes(x = .iteration, y = pr_param, colour = factor(.chain)),
            run_data_orig) +
  geom_step(aes(x = .iteration - 1, y = pr_param),
            run_data) +
  
  coord_cartesian(xlim = c(0, 1000))


ggplot() +
  
  geom_step(aes(x = .iteration, y = eps, colour = factor(.chain)),
            run_data_orig) +
  geom_line(aes(x = .iteration - 1, y = eps, colour = "My code"),
            run_data) +
  
  scale_y_log10() +
  
  coord_cartesian(xlim = c(0, 500))

ggplot() +
  geom_line(aes(x = .iteration, y = eps, colour = "My code"),
            run_data) +
  
  geom_line(aes(x = .iteration, y = eps, colour = "Original code"),
            tibble(.iteration = (1:length(e$accept_rate_tab)) / 2, eps = e$eps_tab)) +
  
  scale_y_log10()


library(targets)
tar_source()

source("replication_paper/common.R")



ix_chain <- 1

chains <- map(
  1:8,
  function(ix_chain) {
    load(
      str_c("../flu-model-kucharski/sero_model/posterior_sero_runs/outputR_f2007_2008_2009_2010_2011_2012_s", 
            ix_chain, "_H3HN_linTRUE.RData")
      ,
      envir = e
    )
    
    
    e$thetatab |> 
      `colnames<-`(c("mu_long", "tau", "tau_2", "omega", "sigma_long", "mu_short", "obs_sd", "disp_k", "sigma_short")) |> 
      as_tibble() |> 
      mutate(.iteration = row_number(),
             .chain = ix_chain,
             .draw = (.iteration - min(.iteration)) + (.chain - 1) * (max(.iteration) - min(.iteration) + 1) + 1,
             .before = 1)
  }
) |> 
  bind_rows()

plot_data <- chains |>
  
  pivot_longer(-c(.iteration, .chain, .draw),
               names_to = "variable")


ggplot() +
  geom_line(aes(x = .iteration, y = value, colour = factor(.chain)),
            plot_data |> filter(.iteration < 10000, variable == "mu_long" | variable == "mu_short")) +
  
  facet_wrap(~variable, scales = "free_y", ncol = 1) +
  
  plot_theme_paper

chains |> 
  filter(.iteration == 1) |> 
  View()

