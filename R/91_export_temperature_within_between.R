# Export supplementary model results and prediction intervals from saved fits.
# Fits and individual records remain private; only aggregate outputs are saved.
.libPaths(c(".Rlib", .libPaths()))
suppressPackageStartupMessages({ library(glmmTMB); library(MASS) })
set.seed(20261026)
base <- "models/within_between_20261026"
out <- paste0("tables/within_between_20261026", Sys.getenv("WB_OUTPUT_SUFFIX"))
dir.create(out, recursive = TRUE, showWarnings = FALSE)
put <- function(x, name) write.csv(x, file.path(out, name), row.names = FALSE)
all_curves <- all_effects <- all_checks <- random_sd <- list()
traits <- strsplit(Sys.getenv("WB_TRAITS", "clutch,fledging"), ",", fixed = TRUE)[[1]]
models <- file.path(base, paste0(traits, "_within_between_selected.rds"))
before <- tools::md5sum(models)

cmp_pmf <- function(mu, phi) {
  nu <- 1 / phi; k <- 0:max(40, ceiling(max(mu) + 20)); lf <- lgamma(k + 1)
  prob <- function(a) {
    z <- outer(a, k) - outer(nu, lf)
    z <- exp(z - do.call(pmax, as.data.frame(z)))
    z / rowSums(z)
  }
  lo <- nu * log(pmax(mu - 3, .001)) - 20
  hi <- nu * log(mu + 3) + 20
  for (i in 1:65) {
    mid <- (lo + hi) / 2; p <- prob(mid); low <- drop(p %*% k) < mu
    lo[low] <- mid[low]; hi[!low] <- mid[!low]
  }
  p <- prob((lo + hi) / 2)
  stopifnot(max(abs(drop(p %*% k) - mu)) < 1e-7,
            max(p[, ncol(p)]) < 1e-10,
            max(abs(drop(p %*% k^2) - mu^2 - compois()$variance(mu, phi))) < 1e-5)
  p
}

for (trait in traits) {
  set.seed(20261026 + match(trait, c("clutch", "fledging")))
  cat(trait, "prediction export\n"); flush.console()
  m <- readRDS(file.path(base, paste0(trait, "_within_between_selected.rds")))
  d <- readRDS(file.path(base, paste0(trait, "_analysis_data.rds")))
  stopifnot(m$fit$convergence == 0, isTRUE(m$sdr$pdHess), nobs(m) == nrow(d))
  b <- fixef(m)
  terms <- c(names(b$cond), paste0("disp~", names(b$disp)),
             if (length(b$zi)) paste0("zi~", names(b$zi)))
  V <- vcov(m, full = TRUE)[terms, terms, drop = FALSE]
  point <- c(unname(b$cond), unname(b$disp), unname(b$zi))
  stopifnot(min(eigen(V, symmetric = TRUE, only.values = TRUE)$values) > 0)
  draws <- MASS::mvrnorm(4000, point, V)
  nc <- length(b$cond); nd <- length(b$disp)
  design <- if (trait == "clutch")
    ~ rel_LD + temp_within + temp_between + year_c + age_class else
    ~ rel_LD + age_class + temp_within + temp_between + year_c
  centre <- unique(round(d$temp_mean - d$temp_mean_c, 10))
  reference <- function(n) data.frame(
    rel_LD = 0, temp_within = rep(0, n), temp_between = 0, year_c = 0,
    age_class = factor(rep("YOUNG", n), levels = c("YOUNG", "OLD")))
  evaluate <- function(newdata) {
    X <- model.matrix(design, newdata)
    location <- drop(X %*% b$cond); scale <- drop(X %*% b$disp)
    locations <- X %*% t(draws[, seq_len(nc), drop = FALSE])
    scales <- X %*% t(draws[, nc + seq_len(nd), drop = FALSE])
    phi <- exp(scale); phis <- exp(scales)
    if (trait == "clutch") {
      mu <- exp(location); mus <- exp(locations)
      consistency <- sqrt(compois()$variance(mu, phi))
      consistencies <- matrix(sqrt(compois()$variance(as.vector(mus), as.vector(phis))),
                               nrow(newdata), nrow(draws))
    } else {
      zero <- plogis(b$zi[1]); zeros <- plogis(draws[, nc + nd + 1])
      mu <- (1 - zero) * plogis(location)
      mus <- sweep(plogis(locations), 2, 1 - zeros, "*")
      consistency <- phi; consistencies <- phis
    }
    native_data <- newdata; native_data$female_id <- NA; native_data$year_f <- NA
    native_mean <- predict(m, newdata = native_data, type = "response", re.form = NA,
                           allow.new.levels = TRUE)
    native_phi <- predict(m, newdata = native_data, type = "disp", re.form = NA,
                          allow.new.levels = TRUE)
    stopifnot(max(abs(native_mean - mu)) < 1e-7,
              max(abs(native_phi - phi)) < 1e-7)
    list(mean = mu, means = mus, consistency = consistency, consistencies = consistencies)
  }
  for (predictor in c("temp_within", "temp_between")) {
    range <- quantile(d[[predictor]], c(.05, .95))
    newdata <- reference(101)
    newdata[[predictor]] <- seq(range[1], range[2], length.out = 101)
    pred <- evaluate(newdata)
    for (component in c("mean", "consistency")) {
      values <- pred[[if (component == "mean") "means" else "consistencies"]]
      ci <- t(apply(values, 1, quantile, c(.025, .975)))
      all_curves[[paste(trait, predictor, component)]] <- data.frame(
        trait = trait, predictor = predictor, component = component,
        x = newdata[[predictor]] + if (predictor == "temp_between") centre else 0,
        estimate = pred[[component]], lower = ci[, 1], upper = ci[, 2])
    }
    contrast <- reference(2); contrast[[predictor]] <- c(0, 1)
    z <- evaluate(contrast)
    mean_change <- z$means[2, ] - z$means[1, ]
    ratio <- z$consistencies[2, ] / z$consistencies[1, ]
    all_effects[[paste(trait, predictor)]] <- data.frame(
      trait = trait, predictor = predictor,
      mean_change_per_C = z$mean[2] - z$mean[1],
      mean_lower = quantile(mean_change, .025), mean_upper = quantile(mean_change, .975),
      consistency_measure = if (trait == "clutch") "conditional_SD" else "component_precision",
      consistency_ratio_per_C = z$consistency[2] / z$consistency[1],
      ratio_lower = quantile(ratio, .025), ratio_upper = quantile(ratio, .975))
  }
  re <- VarCorr(m)
  random_sd[[trait]] <- do.call(rbind, lapply(names(re), function(component)
    do.call(rbind, lapply(names(re[[component]]), function(group)
      data.frame(trait = trait, component = component, group = group,
                 SD = unname(attr(re[[component]][[group]], "stddev")))))))

  # Conditional plug-in checks at fitted group modes; no parameter uncertainty
  # or held-out assessment is implied by these simulation intervals.
  mu <- predict(m, type = if (trait == "clutch") "response" else "conditional")
  phi <- predict(m, type = "disp")
  statistic <- if (trait == "clutch") function(y)
    c(mean = mean(y), variance = var(y), fraction_six = mean(y == 6),
      outside_2_9 = mean(y < 2 | y > 9)) else function(y)
    c(mean_fraction = mean(y / d$clutch), variance_fraction = var(y / d$clutch),
      complete_failure = mean(y == 0), complete_success = mean(y == d$clutch))
  observed <- statistic(if (trait == "clutch") d$clutch else d$fledged)
  if (trait == "clutch") {
    pmf <- cmp_pmf(mu, phi)
    cdf <- t(apply(pmf, 1, cumsum))
    simulation <- replicate(300, statistic(rowSums(cdf < runif(nrow(d)))))
  } else {
    zero <- predict(m, type = "zprob")
    simulation <- replicate(300, {
      y <- rbinom(nrow(d), d$clutch, rbeta(nrow(d), mu * phi, (1 - mu) * phi))
      y[runif(nrow(d)) < zero] <- 0
      statistic(y)
    })
  }
  interval <- t(apply(simulation, 1, quantile, c(.025, .975)))
  all_checks[[trait]] <- data.frame(trait = trait, statistic = names(observed),
                                    observed = unname(observed), predicted = rowMeans(simulation),
                                    lower = interval[, 1], upper = interval[, 2], simulations = 300)
  put(do.call(rbind, all_curves), "prediction_curves.csv")
  put(do.call(rbind, all_effects), "reference_effects_per_degree.csv")
  put(do.call(rbind, random_sd), "random_effect_SD.csv")
  put(do.call(rbind, all_checks), "conditional_simulation_checks.csv")
  cat(trait, "export complete\n"); flush.console()
}
stopifnot(identical(before, tools::md5sum(models)))
capture.output(sessionInfo(), file = file.path(out, "export_session_info.txt"))
saveRDS(list(fixed_effects = read.csv(file.path(out, "fixed_effects.csv")),
             contrasts = read.csv(file.path(out, "within_between_contrasts.csv")),
             curves = do.call(rbind, all_curves),
             reference_effects = do.call(rbind, all_effects)),
        file.path(out, "aggregate_results.rds"), version = 3)
