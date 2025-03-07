
library(tidyverse)
library(arrow)

library(tidybayes)


read_draws <- function(filename) {
  read_parquet(filename) %>%
    rename(.chain = chain, .draw = iteration) %>% 
    mutate(.iteration = 1)
}

