
library(tidyverse)
library(arrow)

library(tidybayes)
library(bayesplot)


read_draws <- function(filename) {
  read_parquet(filename) %>%
    rename(.chain = chain, .iteration = iteration) %>% 
    mutate(Chain = .chain, # Copy for whatever reason.
           .draw = (.iteration - min(.iteration)) + (.chain - 1) * (max(.iteration) - min(.iteration) + 1),
           .before = 3) %>%
    rename_with(function(x) str_remove(x, " ")) 
}

