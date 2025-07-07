

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

reformat_pigeons_chain <- function(chain_pigeons, model_data) {
  n_subjects <- length(model_data$age_distribution)
  n_t_steps <- length(model_data$modelled_years)
  
  inf_columns_new <- colnames(chain_pigeons) |> 
    keep(~ str_detect(.x, "infections")) |> 
    map_dbl(~ as.numeric(str_extract(.x, "\\d+")) - 1) |> 
    map_chr(~ str_c("infections[", .x %% n_t_steps + 1, ",", .x %/% n_t_steps + 1, "]"))
  
  col_names <- colnames(chain_pigeons)
  col_names[str_detect(col_names, "infections")] <- inf_columns_new
  
  colnames(chain_pigeons) <- col_names
  
  return(chain_pigeons)
}

get_lp_mixis <- function(
    chain, n_warmup, model_data,
    turing_model_name, infection_prior, fixed_params,
    mixture_importance_sampling = TRUE, add_name = NULL
) {
  if(!mixture_importance_sampling) {
    return(tibble(lp_mixis = NA, name = add_name, .chain = NA))
  }
  
  chain <- chain %>%
    filter(iteration > n_warmup)
  
  logp <- get_julia_function("pointwise_likelihood")(
    chain, model_data, turing_model_name, infection_prior,
    fixed_params = fixed_params
  )
  
  require(matrixStats)
  
  map(
    1:max(chain$chain),
    function(ix_chain) {
      rows_chain <- chain$chain == ix_chain
      
      l_common_mix <- rowLogSumExps(-logp[rows_chain, ])
      log_weights <- -logp[rows_chain, ] - l_common_mix
      elpd_mixis <- logSumExp(-l_common_mix) - rowLogSumExps(t(log_weights))
      
      tibble(lp_mixis = sum(elpd_mixis), .chain = ix_chain, name = add_name)
    }
  ) %>%
    bind_rows()
}


