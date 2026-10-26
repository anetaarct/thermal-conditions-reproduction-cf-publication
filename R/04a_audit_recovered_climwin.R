# Recovery audit, added 26 October 2026. This is separate from the unchanged
# historical script R/04_climwin_laying_date.R. It does not repeat a window
# search or randomisation. It reads aggregate caches and fits three small LMs.
# Run from the publication repository root with base R (no extra packages).

cache <- readRDS("cache/climwin/recovered_search_data.rds")
annual <- read.csv("data_derived/annual_raw_median_laying_date.csv")
daily <- read.csv("data_derived/smhi_hoburg_daily_temperature.csv")
daily$date <- as.Date(daily$date)
stopifnot(nrow(annual) == 45L, !anyDuplicated(daily$date),
          all(diff(daily$date) == 1), nrow(cache$randomised) == 500L)
out <- "tables/climwin_recovered"
dir.create(out, recursive = TRUE, showWarnings = FALSE)
put <- function(x, name) write.csv(x, file.path(out, name), row.names = FALSE)

best_row <- function(x) {
  x <- x[x$WindowOpen - x$WindowClose >= 2, , drop = FALSE]
  x[which.min(x$deltaAICc), , drop = FALSE]
}
best <- lapply(cache[c("main", "first23", "second22")],
               function(x) best_row(x$Dataset))

# Recover complete-week exposure without installing climwin. In these cached
# searches the seven-day bins start on 1 January of each calendar year.
daily_year <- as.integer(format(daily$date, "%Y"))
daily_week <- (as.integer(format(daily$date, "%j")) - 1L) %/% 7L
exposure <- function(y, open, close) {
  ref <- as.Date(sprintf("%d-05-07", y))
  ref_week <- (as.integer(format(ref, "%j")) - 1L) %/% 7L
  keep <- daily_year == y & daily_week >= ref_week - open &
    daily_week <= ref_week - close
  stopifnot(sum(keep) == (open - close + 1L) * 7L)
  data.frame(year = y, climate = mean(daily$temperature[keep]),
             first_date = min(daily$date[keep]), last_date = max(daily$date[keep]),
             n_days = sum(keep))
}
window_exposure <- function(b, years = annual$year) {
  do.call(rbind, lapply(years, exposure,
                       open = b$WindowOpen, close = b$WindowClose))
}

selection <- do.call(rbind, lapply(names(best), function(label) {
  b <- best[[label]]
  ordinary <- window_exposure(b, 1981L)
  leap <- window_exposure(b, 1980L)
  data.frame(analysis = label, open_weeks = b$WindowOpen,
             close_weeks = b$WindowClose, delta_AICc = b$deltaAICc,
             temperature_beta = b$ModelBeta, SE = b$Std.Error,
             non_leap_start = format(ordinary$first_date, "%d %b"),
             non_leap_end = format(ordinary$last_date, "%d %b"),
             leap_start = format(leap$first_date, "%d %b"),
             leap_end = format(leap$last_date, "%d %b"),
             n_days = ordinary$n_days)
}))
put(selection, "selected_windows_audit.csv")
put(cache$main$Dataset, "candidate_windows.csv")
put(cache$randomised, "randomised_windows.csv")

null <- cache$randomised$deltaAICc
observed <- best$main$deltaAICc
exceed <- sum(null <= observed)
randomisation <- data.frame(
  observed_delta_AICc = observed, repeats = length(null),
  exceedances = exceed, plus_one_p = (exceed + 1) / (length(null) + 1),
  null_min = min(null), null_max = max(null))
put(randomisation, "randomisation_audit.csv")
stopifnot(exceed == 0L, abs(observed + 38.60159858397057) < 1e-8)

main_exposure <- window_exposure(best$main)
training_window <- window_exposure(best$first23)
second_exposure <- window_exposure(best$second22)
annual$main_climate <- main_exposure$climate
annual$climate <- training_window$climate
annual$week_start <- training_window$first_date
annual$week_end <- training_window$last_date
annual$week_n_days <- training_window$n_days
annual$partition <- ifelse(annual$year <= 2002, "training", "held_out")
put(annual, "annual_weekly_exposure.csv")

main_fit <- lm(median_LD ~ main_climate, data = annual)
train <- subset(annual, year <= 2002)
test <- subset(annual, year >= 2003)
training_fit <- lm(median_LD ~ climate, data = train)
held_out_fit <- lm(median_LD ~ climate, data = test)
second <- subset(annual, year >= 2003)
second$climate <- second_exposure$climate[annual$year >= 2003]
second_fit <- lm(median_LD ~ climate, data = second)

# These agreements verify the bin reconstruction against the original search
# output, rather than assuming that a similar-looking daily window is enough.
stopifnot(abs(coef(main_fit)[2] - best$main$ModelBeta) < 1e-10,
          abs(coef(training_fit)[2] - best$first23$ModelBeta) < 1e-10,
          abs(coef(second_fit)[2] - best$second22$ModelBeta) < 1e-10)
aicc <- function(m) {
  k <- attr(logLik(m), "df"); n <- nobs(m)
  AIC(m) + 2 * k * (k + 1) / (n - k - 1)
}
reconstructed_delta <- aicc(main_fit) - aicc(lm(median_LD ~ 1, annual))
stopifnot(abs(reconstructed_delta - observed) < 1e-8)

ci <- confint(held_out_fit)["climate", ]
predicted <- predict(training_fit, newdata = test)
legacy <- cache$legacy_singlewin
stopifnot(isTRUE(all.equal(as.numeric(legacy$model_data$yvar),
                          as.numeric(test$median_LD))))
validation <- data.frame(
  training_seasons = nrow(train), held_out_seasons = nrow(test),
  training_beta = unname(coef(training_fit)[2]),
  held_out_beta = unname(coef(held_out_fit)[2]),
  held_out_CI_low = ci[1], held_out_CI_high = ci[2],
  held_out_R_squared = summary(held_out_fit)$r.squared,
  held_out_p = coef(summary(held_out_fit))["climate", "Pr(>|t|)"],
  prediction_RMSE = sqrt(mean((predicted - test$median_LD)^2)),
  baseline_RMSE = sqrt(mean((mean(train$median_LD) - test$median_LD)^2)),
  prediction_bias = mean(predicted - test$median_LD),
  legacy_singlewin_beta = unname(legacy$coefficients["climate"]),
  max_legacy_exposure_difference = max(abs(legacy$model_data$climate - test$climate)))
rownames(validation) <- NULL
put(validation, "temporal_validation_audit.csv")
put(data.frame(year = test$year, median_LD = test$median_LD,
               complete_week_climate = test$climate,
               legacy_singlewin_climate = legacy$model_data$climate),
    "legacy_exposure_comparison.csv")
stopifnot(abs(validation$held_out_beta + 1.565527) < 1e-6,
          abs(validation$prediction_RMSE - 3.186713) < 1e-6)

saveRDS(list(selection = selection, randomisation = randomisation,
             temporal_validation = validation),
        "cache/climwin/recovery_audit_results.rds", version = 3)
capture.output(sessionInfo(), file = file.path(out, "recovery_session_info.txt"))
print(selection, row.names = FALSE)
print(randomisation, row.names = FALSE)
print(validation, row.names = FALSE)
