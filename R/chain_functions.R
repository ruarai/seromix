
get_inf_data <- function(chain, model_data) {
  chain %>%
    
    spread_draws(infections[ix_t, ix_subject]) %>%
    
    group_by(ix_t, ix_subject) %>%
    summarise(p = 1 - sum(infections == 0) / n()) %>%
    left_join(model_data$infections %>% mutate(inf = TRUE)) %>%
    mutate(inf = replace_na(inf, FALSE)) %>%
    
    left_join(model_data$subject_birth_data) %>%
    mutate(ix_t_birth = replace_na(ix_t_birth, 0)) %>% 
    filter(ix_t >= ix_t_birth)
}


get_inf_accuracy <- function(chain, model_data) {
  chain %>%
    filter(.iteration %% 5 == 0) %>% 
    
    spread_draws(infections[ix_t, ix_subject]) %>%
    
    left_join(model_data$infections %>% mutate(inf = TRUE)) %>%
    mutate(inf = replace_na(inf, FALSE),
           infection = infections > 0) %>%
    
    left_join(model_data$subject_birth_data) %>% 
    filter(ix_t >= ix_t_birth) %>%
    group_by(.iteration, .chain) %>%
    summarise(accuracy = sum(inf == infection) / n())
}
