files <- c(
  clutch_mean = "tables/temporal_models_20261022/clutch_year_mean.csv",
  clutch_sd = "tables/temporal_models_20261022/clutch_year_sd.csv",
  fledging_success = "tables/temporal_models_20261022/fledging_year_success.csv",
  fledging_precision = "tables/temporal_models_20261022/fledging_year_phi.csv"
)
par(mfrow = c(2, 2), mar = c(4, 4, 2, 1))
for (name in names(files)) {
  d <- read.csv(files[[name]])
  plot(d$x, d$estimate, type = "n", xlab = "Calendar year",
       ylab = name, main = gsub("_", " ", name), bty = "l",
       ylim = range(d$CI_low, d$CI_high, finite = TRUE))
  polygon(c(d$x, rev(d$x)), c(d$CI_low, rev(d$CI_high)),
          border = NA, col = adjustcolor("#0072B2", .18))
  lines(d$x, d$estimate, col = "#0072B2", lwd = 2)
}
