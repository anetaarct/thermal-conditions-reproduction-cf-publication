d <- read.csv("data_derived/annual_main_window_temperature.csv")
t <- (d$year - mean(d$year))/10
y <- d$temp_mean
start <- lm(y ~ t)
# The same Gaussian location-scale likelihood, fitted using base R.
nll <- function(b) {
  mu <- b[1] + b[2] * t
  log_sd <- b[3] + b[4] * t
  sum(log_sd + 0.5 * ((y - mu) / exp(log_sd))^2 + log(2*pi)/2)
}
fit <- optim(c(coef(start), log(sqrt(mean(residuals(start)^2))), 0),
             nll, method = "BFGS", hessian = TRUE,
             control = list(reltol = 1e-12, maxit = 2000))
stopifnot(fit$convergence == 0, min(eigen(fit$hessian, symmetric = TRUE)$values) > 0)
se <- sqrt(diag(solve(fit$hessian)))
print(data.frame(term = c("mean_intercept", "warming_per_decade", "log_SD_intercept", "log_SD_per_decade"),
  estimate = fit$par, lower = fit$par - qnorm(.975)*se, upper = fit$par + qnorm(.975)*se), row.names = FALSE)
