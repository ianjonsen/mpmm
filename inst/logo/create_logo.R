require(hexSticker, quietly = TRUE)
require(mpmm)

d <- readRDS("inst/logo/data.RDS")
fit <- mpmm(~ ice + chl + (ice + chl | id), data = d, mpmm_control(REML = TRUE))

p <- plot(fit,lwd = c(0.1, 0.25)) + theme_minimal() +
  labs(title = expression(gamma[t]*" ~ ice + chl + (ice + chl | id)")) +
 #   "\u03b3", " ~ ice + chl + (ice + chl | id)")) +
  theme(strip.text.x = element_blank(),
        axis.title.y = element_blank(),
        panel.grid = element_line(size=0.05, color = "grey70"),
        plot.title = element_text(hjust=0.45, vjust=0, size = 16, color = "#045a8d"))

s <- sticker(
  p,
  package = "mpmm",
  p_size = 24,
  p_y = 1.6,
  p_family = "sans",
  #p_fontface = "bold",
  p_color = "#045a8d",
  h_color = "#045a8d",
  s_x = 0.8,
  s_y = 0.8,
  s_width = 1.65,
  s_height = 1.35,
  h_fill =  "white",
  spotlight = FALSE,
  l_x = 0.94,
  l_y = 1.08,
  l_width = 2,
  filename = "inst/logo/mpmm_logo.png"
)

