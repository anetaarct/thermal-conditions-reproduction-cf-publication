curves <- read.csv("tables/joint_calendar_year_20261027/curves.csv")
par(mfrow = c(2,2), mar = c(4,4,2,1))
colours <- c(baseline = "#0072B2", within_between = "#009E73")
for (trait in c("clutch", "fledging")) for (component in c("mean", "consistency")) {
  d <- subset(curves, curves$trait == trait & curves$component == component)
  ylabel <- if (component == "mean") {
    if (trait == "clutch") "Expected clutch size" else "Expected fledging success"
  } else if (trait == "clutch") "Conditional clutch-size SD" else "Component precision"
  plot(range(d$year), range(d$lower, d$upper), type = "n", xlab = "Calendar year", ylab = ylabel,
       main = paste(trait, component), bty = "l")
  for (model in names(colours)) {
    z <- subset(d, d$model == model)
    polygon(c(z$year, rev(z$year)), c(z$lower, rev(z$upper)), border = NA,
            col = adjustcolor(colours[model], alpha.f = .13))
    lines(z$year, z$estimate, col = colours[model], lwd = 2,
          lty = if (model == "baseline") 1 else 2)
  }
  legend("topleft", c("Single temperature", "Within-between"), col = colours,
         lty = c(1,2), lwd = 2, bty = "n", cex = .7)
}
