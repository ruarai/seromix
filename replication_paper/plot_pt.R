

chains_pt <- bind_rows(
  arrow::read_parquet("runs/hanam_2018/chain_pigeons_6_1.parquet") |> mutate(chain = 1),
  arrow::read_parquet("runs/hanam_2018/chain_pigeons_6_2.parquet") |> mutate(chain = 2),
) |> 
  reformat_pigeons_chain(tar_read(hanam_2018)) |> 
  make_chain_subset(tar_read(hanam_2018), NULL)



ggplot(chains_pt) +
  geom_line(aes(x = .iteration, y = total_inf, colour = factor(.chain)))+
  
  scale_colour_discrete(type = rep(RColorBrewer::brewer.pal(9, "Blues")[c(6,9)], 10)) +
  
  
  plot_theme_paper +
  
  theme(strip.text = element_markdown(size = 12, family = "Utopia"),
        panel.grid.major.y = element_gridline,
        legend.position = "none") +
  
  ggtitle("Ha Nam study inference results")
