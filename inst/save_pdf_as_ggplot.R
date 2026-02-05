library(magick)
img <- image_read("figures/model_schematic_simple.pdf")
grob=rasterGrob(image_draw(img))

p=ggplot() +
  # Add other ggplot layers as needed (e.g., geom_point, geom_line)
  # For a blank plot, you can use geom_blank()
  geom_blank() +
  annotation_custom(grob, xmin = -Inf, xmax = Inf, ymin = -Inf, ymax = Inf) +
  # Adjust xlim and ylim if you want to place the image within specific coordinates
  # For example: xlim(0, 10) + ylim(0, 5)
  theme_void() # Or other theme options

save(p,file="figures/model_schematic_simple.RData")
plot_grid(p)
