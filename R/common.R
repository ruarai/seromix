
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



save_hdf5 <- function(data_list, filename) {
  if(file.exists(filename)) {
    file.remove(filename)
  }
  
  for (i in 1:length(data_list)) {
    rhdf5::h5write(data_list[[i]], filename, names(data_list)[[i]])
  }
}
