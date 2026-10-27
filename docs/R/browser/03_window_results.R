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
