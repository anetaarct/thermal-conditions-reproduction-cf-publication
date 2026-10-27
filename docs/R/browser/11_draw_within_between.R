curves <- read.csv("tables/within_between_20261026/prediction_curves.csv")
trait <- "clutch" # Change to "fledging" and run again.
par(mfrow = c(2,2), mar = c(4,4,2,1))
for (component in c("mean", "consistency")) for (predictor in c("temp_within", "temp_between")) {
  d <- subset(curves, curves$trait == trait & curves$component == component & curves$predictor == predictor)
  col <- if (component == "mean") "#0072B2" else "#D55E00"
  plot(d$x, d$estimate, type = "n", ylim = range(d$lower, d$upper),
       xlab = predictor, ylab = component, main = paste(trait, predictor), bty = "l")
  polygon(c(d$x, rev(d$x)), c(d$lower, rev(d$upper)), border = NA, col = adjustcolor(col, .2))
  lines(d$x, d$estimate, col = col, lwd = 2)
}
