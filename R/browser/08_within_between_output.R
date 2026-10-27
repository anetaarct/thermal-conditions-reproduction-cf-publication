results <- readRDS("tables/within_between_20261026/aggregate_results.rds")
print(subset(results$fixed_effects, model == "within_between" &
             term %in% c("temp_within", "temp_between")), row.names = FALSE)
print(results$contrasts, row.names = FALSE)
effects <- results$reference_effects
effects$consistency_change_percent <- 100 * (effects$consistency_ratio_per_C - 1)
print(effects[, c("trait", "predictor", "mean_change_per_C", "consistency_measure",
                  "consistency_change_percent", "ratio_lower", "ratio_upper")], row.names = FALSE)
