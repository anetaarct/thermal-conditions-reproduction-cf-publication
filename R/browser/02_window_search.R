annual <- read.csv("data_derived/annual_raw_median_laying_date.csv")
daily <- read.csv("data_derived/smhi_hoburg_daily_temperature.csv")
saved <- read.csv("tables/climwin_recovered/selected_windows_audit.csv")
daily$date <- as.Date(daily$date)
daily$week <- (as.integer(format(daily$date, "%j")) - 1L) %/% 7L
windows <- subset(expand.grid(open = 0:17, close = 0:17), open - close >= 2)
exposure <- matrix(NA_real_, nrow(annual), nrow(windows))
for (i in seq_len(nrow(annual))) {
  y <- annual$year[i]
  w <- subset(daily, year == y)
  ref <- (as.integer(format(as.Date(sprintf("%d-05-07", y)), "%j")) - 1L) %/% 7L
  for (j in seq_len(nrow(windows))) {
    keep <- w$week >= ref - windows$open[j] & w$week <= ref - windows$close[j]
    stopifnot(sum(keep) == 7 * (windows$open[j] - windows$close[j] + 1))
    exposure[i, j] <- mean(w$temperature[keep])
  }
}
# Gaussian linear-model AICc: intercept+slope+residual SD = 3 parameters;
# the intercept-only comparison has 2 parameters.
n <- nrow(annual)
X <- scale(exposure)
y <- as.numeric(scale(annual$median_LD))
r <- drop(crossprod(X, y)) / (n - 1)
delta <- n * log1p(-pmin(r^2, 1 - 1e-15)) + 2 + 24/(n-4) - 12/(n-3)
best <- which.min(delta)
fit <- lm(annual$median_LD ~ exposure[, best])
answer <- data.frame(open_weeks = windows$open[best], close_weeks = windows$close[best],
                     delta_AICc = delta[best], temperature_beta = coef(fit)[2])
stopifnot(answer$open_weeks == saved$open_weeks[1],
          answer$close_weeks == saved$close_weeks[1],
          abs(answer$delta_AICc - saved$delta_AICc[1]) < 1e-5,
          abs(answer$temperature_beta - saved$temperature_beta[1]) < 1e-7)
print(answer, row.names = FALSE)
# Re-run complete eligible-window searches under randomisation.
# This is a new reproducible set of shuffles, not a replay of saved shuffles.
set.seed(20261027)
B <- 500L
permuted <- replicate(B, sample(y))
correlations <- crossprod(X, permuted) / (n - 1)
null <- n * log1p(-apply(correlations^2, 2, max)) + 2 + 24/(n-4) - 12/(n-3)
cat("Eligible windows:", nrow(windows), "\nRandomised full searches:", B,
    "\nExceedances:", sum(null <= delta[best]),
    "\nPlus-one p:", (1 + sum(null <= delta[best]))/(B + 1), "\n")
