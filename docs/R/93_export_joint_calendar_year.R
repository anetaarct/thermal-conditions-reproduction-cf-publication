# Export continuous-year associations from the two temperature formulations.
# Run from the repository root; fitted objects and records stay in models/.
.libPaths(c(".Rlib", .libPaths()))
suppressPackageStartupMessages({library(glmmTMB); library(MASS)})
base <- "models/within_between_20261026"
out <- "tables/joint_calendar_year_20261027"
dir.create(out, recursive = TRUE, showWarnings = FALSE)
curves <- contrasts <- coefficients <- list()
for (trait in c("clutch", "fledging")) for (model in c("baseline", "within_between")) {
  cat(trait, model, "year export\n"); flush.console()
  set.seed(20261027 + match(trait, c("clutch", "fledging")) + 10L * (model == "within_between"))
  m <- readRDS(file.path(base, paste0(trait, "_", model, "_selected.rds")))
  d <- readRDS(file.path(base, paste0(trait, "_analysis_data.rds")))
  stopifnot(m$fit$convergence == 0, isTRUE(m$sdr$pdHess), nobs(m) == nrow(d))
  b <- fixef(m)
  terms <- c(names(b$cond), paste0("disp~", names(b$disp)),
             if (length(b$zi)) paste0("zi~", names(b$zi)))
  V <- vcov(m, full = TRUE)[terms, terms, drop = FALSE]
  point <- c(unname(b$cond), unname(b$disp), unname(b$zi))
  draws <- MASS::mvrnorm(4000, point, V)
  nc <- length(b$cond); nd <- length(b$disp)
  year_center <- unique(round(d$year - d$year_c, 8))
  stopifnot(length(year_center) == 1)
  years <- seq(min(d$year), max(d$year), length.out = 101)
  ndat <- data.frame(year_c = years - year_center, rel_LD = 0,
    temp_mean_c = 0, temp_within = 0, temp_between = 0,
    age_class = factor("YOUNG", levels = c("YOUNG", "OLD")))
  form <- if (trait == "clutch") {
    if (model == "baseline") ~rel_LD + temp_mean_c + year_c + age_class else
      ~rel_LD + temp_within + temp_between + year_c + age_class
  } else {
    if (model == "baseline") ~rel_LD + age_class + temp_mean_c + year_c else
      ~rel_LD + age_class + temp_within + temp_between + year_c
  }
  X <- model.matrix(form, ndat)
  location <- drop(X %*% b$cond); scale <- drop(X %*% b$disp)
  locations <- X %*% t(draws[, seq_len(nc), drop = FALSE])
  scales <- X %*% t(draws[, nc + seq_len(nd), drop = FALSE])
  phi <- exp(scale); phis <- exp(scales)
  if (trait == "clutch") {
    mean <- exp(location); means <- exp(locations)
    consistency <- sqrt(compois()$variance(mean, phi))
    consistencies <- matrix(sqrt(compois()$variance(as.vector(means), as.vector(phis))), 101, 4000)
  } else {
    mean <- (1 - plogis(b$zi[1])) * plogis(location)
    means <- sweep(plogis(locations), 2, 1 - plogis(draws[, nc + nd + 1]), "*")
    consistency <- phi; consistencies <- phis
  }
  check_data <- ndat[c(1, 101), ]; check_data$female_id <- NA; check_data$year_f <- NA
  stopifnot(max(abs(predict(m, newdata = check_data, type = "response", re.form = NA,
                            allow.new.levels = TRUE) - mean[c(1, 101)])) < 1e-7)
  for (component in c("mean", "consistency")) {
    values <- if (component == "mean") means else consistencies
    ci <- t(apply(values, 1, quantile, c(.025, .975)))
    curves[[paste(trait, model, component)]] <- data.frame(trait, model, component,
      year = years, estimate = if (component == "mean") mean else consistency,
      lower = ci[, 1], upper = ci[, 2])
  }
  difference <- means[101, ] - means[1, ]; ratio <- consistencies[101, ] / consistencies[1, ]
  contrasts[[paste(trait, model)]] <- data.frame(trait, model, first_year = min(years), last_year = max(years),
    mean_first = mean[1], mean_last = mean[101], mean_change = mean[101] - mean[1],
    mean_lower = quantile(difference, .025), mean_upper = quantile(difference, .975),
    consistency_ratio = consistency[101] / consistency[1],
    ratio_lower = quantile(ratio, .025), ratio_upper = quantile(ratio, .975))
  for (component in c("cond", "disp")) {
    z <- summary(m)$coefficients[[component]]["year_c", ]
    coefficients[[paste(trait, model, component)]] <- data.frame(trait, model, component,
      estimate = z[1], SE = z[2], lower = z[1] - qnorm(.975) * z[2], upper = z[1] + qnorm(.975) * z[2])
  }
  write.csv(do.call(rbind, curves), file.path(out, "curves.csv"), row.names = FALSE)
  write.csv(do.call(rbind, contrasts), file.path(out, "contrasts.csv"), row.names = FALSE)
  write.csv(do.call(rbind, coefficients), file.path(out, "year_coefficients.csv"), row.names = FALSE)
}
capture.output(sessionInfo(), file = file.path(out, "session_info.txt"))
cat("Joint calendar-year export complete; no models refitted.\n")
