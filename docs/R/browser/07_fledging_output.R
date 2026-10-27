results <- readRDS("tables/within_between_20261026/aggregate_results.rds")
f <- subset(results$fixed_effects, trait == "fledging" & model == "baseline")
print(f, row.names = FALSE)
b <- subset(f, component == "cond")
b <- setNames(b$estimate, b$term)
zi <- subset(f, component == "zi" & term == "(Intercept)")$estimate
# Reference: yearling, seasonal median date, centred year, random effects zero.
temperature_deviation <- c(0, 1)
eta <- b["(Intercept)"] + b["temp_mean_c"] * temperature_deviation
expected <- (1 - plogis(zi)) * plogis(eta)
print(data.frame(temperature_deviation, expected), row.names = FALSE)
