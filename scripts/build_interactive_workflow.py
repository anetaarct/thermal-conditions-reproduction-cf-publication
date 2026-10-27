"""Build browser-executable base-R steps from public analysis inputs only."""
from pathlib import Path
import html, json, re

root = Path(__file__).resolve().parent.parent
scripts = root / 'R/browser'
scripts.mkdir(parents=True, exist_ok=True)
labs = {}

def lab(key, title, code, note='Runs from public input files.'):
    code = code.strip() + '\n'
    (scripts / (key + '.R')).write_bytes(code.encode('utf-8'))
    assets = sorted(set(re.findall(r'(?:read\.csv|readRDS)\("([^"]+)"', code)) |
                        set(re.findall(r'"((?:data_derived|tables|cache)/[a-zA-Z0-9_./-]+)"', code)))
    for asset in assets:
        assert (root / asset).is_file(), asset
    labs[key] = dict(title=title, note=note, assets=assets, script=f'R/browser/{key}.R')
    block = f'''```{{=html}}
<section class="r-lab" data-r-lab="{key}" aria-label="{html.escape(title)}">
<div class="r-lab-heading"><span class="r-lab-badge">Run R</span><h3>{html.escape(title)}</h3></div>
<p class="r-lab-note">{html.escape(note)}</p>
<label for="r-code-{key}">Editable R code</label>
<textarea id="r-code-{key}" class="r-lab-code" spellcheck="false" rows="{min(22,max(7,len(code.splitlines())))}">{html.escape(code)}</textarea>
<div class="r-lab-actions"><button class="r-run" type="button">Run code</button><button class="r-reset" type="button">Reset code</button><button class="r-download" type="button">Download R code</button><span class="r-status" role="status" aria-live="polite">Ready to run</span></div>
<div class="r-output-label">Output appears here after running.</div><pre class="r-output" tabindex="0"></pre><div class="r-plots"></div>
<noscript>JavaScript is required to run R here. Download the R script and run it locally.</noscript>
</section>
```
'''
    (root / '_includes' / f'lab-{key}.qmd').write_bytes(block.encode('utf-8'))

lab('01_weather_inputs', 'Load and check the weather inputs', r'''
daily <- read.csv("data_derived/smhi_hoburg_daily_temperature.csv")
annual <- read.csv("data_derived/annual_main_window_temperature.csv")
laying <- read.csv("data_derived/annual_raw_median_laying_date.csv")
daily$date <- as.Date(daily$date)
stopifnot(nrow(daily) == 16802, !anyDuplicated(daily$date),
          nrow(annual) == 46, nrow(laying) == 45)
keep <- format(daily$date, "%m-%d") >= "03-26" &
        format(daily$date, "%m-%d") <= "05-07"
x <- daily[keep, ]
rebuilt <- aggregate(temperature ~ year, x, mean)
days <- table(x$year)
stopifnot(all(days == 43),
  max(abs(rebuilt$temperature - annual$temp_mean[match(rebuilt$year, annual$year)])) < 1e-10)
print(head(annual))
cat("Daily records:", nrow(daily), "\nYears:", nrow(rebuilt),
    "\nDays per window:", unique(as.integer(days)), "\nInput checks passed.\n")
''')

lab('02_window_search', 'Recalculate the sliding-window search', r'''
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
''', 'Fits the eligible linear-window search from the daily data, then runs 500 new randomised searches. Base R; no individual bird records are needed.')

lab('03_window_results', 'Read and check the saved window results', r'''
results <- readRDS("cache/climwin/recovered_search_data.rds")
x <- results$main$Dataset
eligible <- subset(x, WindowOpen - WindowClose >= 2)
best <- eligible[which.min(eligible$deltaAICc), ]
print(best[, c("WindowOpen", "WindowClose", "deltaAICc", "ModelBeta")])
null <- results$randomised$deltaAICc
print(c(randomisations = length(null), exceedances = sum(null <= best$deltaAICc),
        plus_one_p = (1 + sum(null <= best$deltaAICc))/(length(null) + 1)))
validation <- readRDS("cache/climwin/recovery_audit_results.rds")
print(validation$temporal_validation)
''', 'Reads the data-only RDS and recalculates selection and the empirical probability from its 500 saved searches.')

lab('04_environment_model', 'Fit spring-temperature location and scale', r'''
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
''', 'Fits the Gaussian mean and log-SD likelihood to 46 annual temperatures. Coefficients for time are per decade; this is not a reproductive model.')

lab('05_predictor_checks', 'Check temperature-year correlations and inspect VIF', r'''
weather <- read.csv("data_derived/annual_main_window_temperature.csv")
samples <- read.csv("tables/temporal_models_20261022/observed_annual.csv")
for (trait in c("clutch", "fledging")) {
  years <- samples$year[samples$trait == trait]
  a <- weather[weather$year %in% years, ]
  cat(trait, "seasons:", nrow(a), "annual Pearson r:", cor(a$year, a$temp_mean), "\n")
}
# These VIFs were calculated on the full brood-level fixed-design matrices.
# Annual aggregates cannot reconstruct the brood-level age/date design.
vif <- read.csv("tables/predictor_diagnostics/fixed_design_vif.csv")
print(vif, row.names = FALSE)
''', 'Recalculates annual correlations. Displays the saved brood-level VIFs; recalculating those VIFs requires the individual analysis records.')

for trait, key, title in [('clutch','06_clutch_output','Read the clutch model and calculate its mean'),
                           ('fledging','07_fledging_output','Read the fledging model and calculate expected success')]:
    response = 'exp(eta)' if trait == 'clutch' else '(1 - plogis(zi)) * plogis(eta)'
    zi = '' if trait == 'clutch' else 'zi <- subset(f, component == "zi" & term == "(Intercept)")$estimate\n'
    lab(key,title, f'''
results <- readRDS("tables/within_between_20261026/aggregate_results.rds")
f <- subset(results$fixed_effects, trait == "{trait}" & model == "baseline")
print(f, row.names = FALSE)
b <- subset(f, component == "cond")
b <- setNames(b$estimate, b$term)
{zi}# Reference: yearling, seasonal median date, centred year, random effects zero.
temperature_deviation <- c(0, 1)
eta <- b["(Intercept)"] + b["temp_mean_c"] * temperature_deviation
expected <- {response}
print(data.frame(temperature_deviation, expected), row.names = FALSE)
''', 'Reads fitted coefficients and recalculates reference means. It does not refit the mixed model; fitting code is shown below.')

lab('08_within_between_output', 'Recalculate within-between summaries', r'''
results <- readRDS("tables/within_between_20261026/aggregate_results.rds")
print(subset(results$fixed_effects, model == "within_between" &
             term %in% c("temp_within", "temp_between")), row.names = FALSE)
print(results$contrasts, row.names = FALSE)
effects <- results$reference_effects
effects$consistency_change_percent <- 100 * (effects$consistency_ratio_per_C - 1)
print(effects[, c("trait", "predictor", "mean_change_per_C", "consistency_measure",
                  "consistency_change_percent", "ratio_lower", "ratio_upper")], row.names = FALSE)
''', 'Uses real aggregate model outputs. A lower clutch SD means greater consistency; a higher fledging component precision means greater consistency.')

# Year outputs are generated by R/93 before this builder runs.
lab('09_joint_year_output', 'Read year effects from both temperature models', r'''
fixed <- read.csv("tables/joint_calendar_year_20261027/year_coefficients.csv")
change <- read.csv("tables/joint_calendar_year_20261027/contrasts.csv")
print(fixed, row.names = FALSE)
print(change, row.names = FALSE)
''', 'Both models include continuous year and thermal predictors together. These are temperature-adjusted year associations, not year-only trends.')

lab('10_draw_joint_year', 'Draw the calendar-year figure', r'''
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
''', 'Draws both temperature formulations from saved prediction tables. Modify the code and rerun to change the figure.')

lab('11_draw_within_between', 'Draw within-between predictions', r'''
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
''', 'Draws real fitted curves and their saved 95% confidence intervals. No model fitting is repeated.')

lab('12_draw_windows', 'Draw window-search diagnostics', r'''
windows <- read.csv("tables/climwin_recovered/candidate_windows.csv")
randomised <- read.csv("tables/climwin_recovered/randomised_windows.csv")
w <- subset(windows, WindowOpen - WindowClose >= 2)
best <- w[which.min(w$deltaAICc), ]
par(mfrow = c(1,2), mar = c(4,4,2,1))
colours <- hcl.colors(100, "Viridis")
shade <- cut(w$deltaAICc, breaks = 100, labels = FALSE)
plot(w$WindowOpen, w$WindowClose, pch = 15, cex = 1.4, col = colours[shade],
     xlab = "Opening offset (weeks)", ylab = "Closing offset (weeks)", main = "Eligible windows")
points(best$WindowOpen, best$WindowClose, pch = 8, cex = 1.7)
hist(randomised$deltaAICc, xlim = range(best$deltaAICc, randomised$deltaAICc),
     main = "500 saved full-search randomisations", xlab = "Delta AICc", col = "#9CC5B3", border = "white")
abline(v = best$deltaAICc, col = "#D55E00", lwd = 2)
''', 'Draws the saved candidate surface and null distribution; the search itself runs in the sliding-window chapter.')

lab('13_draw_temporal_trends', 'Draw calendar-year trends from saved predictions', r'''
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
''', 'Uses saved model predictions and draws a new figure in your browser. No model is refitted.')

(root/'assets').mkdir(exist_ok=True)
(root/'assets/r-labs.json').write_bytes(json.dumps(labs, indent=2).encode())
print('Built',len(labs),'editable R steps and their includes.')
