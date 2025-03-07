
library(tidyverse)
library(arrow)

library(tidybayes)


clean_df <- function(x) {
  rename(x, .chain = chain, .draw = iteration) %>% 
    mutate(.iteration = 1)
}


ppd_df <- read_parquet("data/ppd_df.parquet") %>%
  clean_df()



plot_data <- ppd_df %>%
  spread_draws(titer_obs[t])


ggplot() +
  geom_point(aes(x = t, y = titer_obs, group = .draw),
             plot_data) 
