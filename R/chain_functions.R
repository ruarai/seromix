

summarise_chain <- function(chain, drop_iterations, model_data, by_chain = TRUE) {
  chain_long <- chain %>%
    clean_chain() %>%
    filter(.iteration > drop_iterations) %>%
    add_total_infections(model_data) %>% 
    select(-starts_with("infections")) %>%
    
    pivot_longer(-c(.iteration, .chain, .draw),
                 names_to = "variable")
  
  if(by_chain) {
    chain_long <- chain_long %>% 
      group_by(.chain, variable)
  } else {
    chain_long <- chain_long %>% 
      group_by(variable)
  }
  
  chain_long %>% 
    
    summarise(
      mean = mean(value),
      median = median(value),
      sd = sd(value),
      q95_lower = quantile(value, 0.05 / 2),
      q95_upper = quantile(value, 1 - 0.05 / 2),
      ess_bulk = posterior::ess_bulk(value),
      ess_bulk = posterior::ess_tail(value)
    )
}

add_total_infections <- function(chain, model_data) {
  n_t_steps <- length(model_data$modelled_years)
  
  valid_inf_cols <- model_data$subject_birth_data %>%
    mutate(ix_t_birth = replace_na(ix_t_birth, 1)) %>%
    rowwise() %>%
    mutate(ix_t = list(ix_t_birth:n_t_steps)) %>% 
    unnest(ix_t) %>%
    arrange(ix_t, ix_subject) %>%
    mutate(col_name = str_c("infections[", ix_t, ",", ix_subject, "]")) %>%
    pull(col_name)
  
  chain %>%
    mutate(total_inf = rowSums(across(starts_with("infections"))))
}