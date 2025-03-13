
library(tidyverse)
library(arrow)

library(tidybayes)


read_draws <- function(filename) {
  # Need to make draw = iteration + chain * ... or something
  read_parquet(filename) %>%
    rename(.chain = chain, .iteration = iteration) %>% 
    mutate(.draw = (.iteration - min(.iteration)) + (.chain - 1) * (max(.iteration) - min(.iteration) + 1),
           .before = 3) %>%
    rename_with(function(x) str_remove(x, " ")) 
}

