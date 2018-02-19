map_tracks <- function(d, xlim = c(0, 150), ylim = c(-70, -44)) {

require(ggplot2)
require(mapproj)

countriesHigh <- NULL
data(countriesHigh, package = "rworldxtra", envir = environment())
wm <- suppressMessages(fortify(countriesHigh))
xl <- xlim
yl <- ylim

p <-
  ggplot() + 
  coord_map(projection = "lambert", parameters = c(-65, -45), xlim = xl, ylim = yl) + 
  xlab("Longitude") + ylab("Latitude")

p <- p +
  geom_point(aes(x = lon, y = lat, colour = d$fitted$g),
             data = d$data,
             size = 0.5) +
  viridis::scale_colour_viridis(name = expression(italic(gamma[t])), begin = 0, end = 1, direction = -1) +
  theme_dark() + 
  theme(axis.title.x = element_blank(),
        axis.text.x = element_blank(),
        axis.ticks.x = element_blank()) +
  theme(axis.title.y = element_blank(),
        axis.text.y = element_blank(),
        axis.ticks.y = element_blank()) +
  theme(legend.title = element_text(hjust = 0.5))

p <- p + 
  geom_polygon(
    data = wm,
    aes_string(x = "long", y = "lat", group = "group"),
    fill = grey(0.95)
  ) 

#if(deparse(substitute(d)) == "fit.pel") {
  p <- p + 
  annotate("text", label = "45 S", x = 3, y = -45, size = 2.5, colour = "white") +
  annotate("text", label = "55 S", x = 3, y = -55, size = 2.5, colour = "white") +
  annotate("text", label = "65 S", x = 3, y = -65, size = 2.5, colour = "white") +
  annotate("text", label = "50 E", x = 50, y = -43, size = 2.5, colour = "white") +
  annotate("text", label = "100 E", x = 100, y = -43, size = 2.5, colour = "white")
# } else {
#   p <- p + 
#     annotate("text", label = "50 S", x = 5, y = -50, size = 2.5, colour = "white") +
#     annotate("text", label = "60 S", x = 5, y = -60, size = 2.5, colour = "white") +
#     annotate("text", label = "50 E", x = 50, y = -45, size = 2.5, colour = "white") +
#     annotate("text", label = "100 E", x = 100, y = -45, size = 2.5, colour = "white")
#   
# }

p
}
