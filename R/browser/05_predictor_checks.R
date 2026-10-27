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
