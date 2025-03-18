

hanam_data_processed$subject_id %>% range()


plots <- map(
  1:13,
  function(i) {
    hanam_data_processed %>%
      filter(subject_id %in% (i * 5):((i + 1) * 5 - 1)) %>% 
      ggplot() +
      geom_point(aes(x = strain_year, y = titre),
                 size = 0.7) +
      
      facet_grid(rows = vars(year_sampled), cols = vars(subject_id)) +
      
      coord_cartesian(ylim = c(0, 9)) +
      
      theme_bw()
  }
)

plot(plots[[13]])


pdf("results/data_2018.pdf", width = 10, height = 6)
for (i in 1:length(plots)){
  plot(plots[[i]])
}
dev.off()


pdf("results/data_2018_long.pdf", width = 10, height = 6)
for(i in 1:max(hanam_data_processed$subject_id)) {
  p <- hanam_data_processed %>%
    filter(subject_id == i) %>% 
    ggplot() +
    geom_point(aes(x = year_sampled, y = titre),
               size = 0.7) +
    
    facet_wrap(~strain_year, ncol = 5) +
    
    coord_cartesian(ylim = c(0, 9)) +
    
    theme_bw()
  
  plot(p)
}
dev.off()




