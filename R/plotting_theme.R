

plot_theme_paper <- list(
  theme_minimal(),
  theme(panel.grid = element_blank(),
        axis.line = element_line(),
        axis.ticks = element_line(),
        text = element_text(family = "Helvetica", colour = "black"),
        line = element_line(linewidth = 0.7),
        axis.text.x = element_text(margin = margin(t = 0.3, unit = "cm"), colour = "black"),
        axis.text.y = element_text(margin = margin(r = 0.3, unit = "cm"), colour = "black"),
        plot.margin = margin(t = 0.2, l = 0.2, b = 0.2, r = 0.2, unit = "cm"),
        axis.ticks.length=unit(-0.1, "cm"),
        axis.title.y = element_markdown(),
        axis.title.x = element_markdown(),
        strip.text.x = element_text(hjust = 0),
        plot.title = element_markdown(),
        plot.subtitle = element_markdown()
  )
)

element_gridline <- element_line(colour = "grey50", linetype = "28", linewidth = 0.3)
element_facet_background <- element_rect(fill = "white", colour = "grey50", linewidth = 0.5)
