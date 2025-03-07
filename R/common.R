
library(tidyverse)
library(arrow)

library(tidybayes)


read_draws <- function(filename) {
  # Need to make draw = iteration + chain * ... or something
  read_parquet(filename) %>%
    rename(.chain = chain, .draw = iteration) %>% 
    mutate(.iteration = .draw + (.chain - 1) * max(.draw),
           .before = 3) %>%
    rename_with(function(x) str_remove(x, " ")) 
}

