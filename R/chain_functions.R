

make_chain_subset <- function(chain, model_data, name) {
  chain |> 
    clean_chain() |>
    add_total_infections(model_data) |>
    select(-starts_with("infections")) |>
    mutate(name = name)
}

summarise_chain <- function(chain, drop_iterations, model_data, by_chain = TRUE, add_name = NULL) {
  chain_long <- chain |>
    clean_chain() |>
    filter(.iteration > drop_iterations) |>
    add_total_infections(model_data) |> 
    select(-starts_with("infections")) |>
    
    pivot_longer(-c(.iteration, .chain, .draw),
                 names_to = "variable")
  
  if(by_chain) {
    chain_long <- chain_long |> 
      group_by(.chain, variable)
  } else {
    chain_long <- chain_long |> 
      group_by(variable)
  }
  
  summ_chain <- chain_long |> 
    
    summarise(
      mean = mean(value),
      median = median(value),
      sd = sd(value),
      q95_lower = quantile(value, 0.05 / 2),
      q95_upper = quantile(value, 1 - 0.05 / 2),
      ess_bulk = posterior::ess_bulk(value),
      ess_bulk = posterior::ess_tail(value)
    )
  
  if(is.null(add_name)) {
    return(summ_chain)
  } else{
    return(summ_chain %>% mutate(name = add_name))
  }
}

add_total_infections <- function(chain, model_data) {
  n_t_steps <- length(model_data$modelled_years)
  
  chain |>
    mutate(total_inf = rowSums(across(starts_with("infections"))))
}
