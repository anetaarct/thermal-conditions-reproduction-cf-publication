# Figures from saved aggregate prediction curves; no model fitting.
curves <- read.csv(paste0("tables/within_between_20261026", Sys.getenv("WB_OUTPUT_SUFFIX"), "/prediction_curves.csv"))
for (trait in unique(curves$trait)) {
  draw <- function() {
    par(mfrow = c(2, 2), mar = c(4.4, 4.4, 2.3, 1), mgp = c(2.7, .8, 0))
    letter <- 0L
    for (component in c("mean", "consistency")) {
      for (predictor in c("temp_within", "temp_between")) {
        letter <- letter + 1L
        d <- curves[curves$trait == trait & curves$component == component &
                      curves$predictor == predictor, ]
        colour <- if (component == "mean") "#0072B2" else "#D55E00"
        xlabel <- if (predictor == "temp_within")
          "Temperature minus female mean (\u00b0C)" else "Female mean temperature (\u00b0C)"
        ylabel <- if (component == "mean") {
          if (trait == "clutch") "Expected clutch size (eggs)" else "Expected egg-to-fledging success"
        } else if (trait == "clutch") "Conditional clutch-size SD" else "Beta-binomial component precision"
        pair <- curves[curves$trait == trait & curves$component == component, ]
        plot(d$x, d$estimate, type = "n", ylim = range(c(pair$lower, pair$upper)),
             xlab = xlabel, ylab = ylabel, bty = "l")
        polygon(c(d$x, rev(d$x)), c(d$lower, rev(d$upper)),
                col = adjustcolor(colour, alpha.f = .18), border = NA)
        lines(d$x, d$estimate, col = colour, lwd = 2)
        title(main = paste(LETTERS[letter], if (predictor == "temp_within") "Within female" else "Between females"),
              adj = 0, cex.main = 1)
      }
    }
  }
  stem <- paste0("figures/Fig_within_between_", trait)
  png(paste0(stem, ".png"), width = 1700, height = 1350, res = 180)
  draw(); dev.off()
  pdf(paste0(stem, ".pdf"), width = 9.44, height = 7.5)
  draw(); dev.off()
}
