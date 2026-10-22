# Read-only diagnostic of the exact data stored in the final Bayesian M1 fit.
# Run from /home/arct/CF_model_comparison_20261013 on Totoro. No model fitting.
model_path <- 'fledging_zi_additive_20261016/M1_zero_inflated_beta_binomial.rds'
before <- unname(tools::md5sum(model_path))
m <- readRDS(model_path)
d <- m$data
stopifnot(nrow(d) == 9790L, length(unique(d$female_id)) == 7132L,
          length(unique(d$year_f)) == 39L,
          all(c('rel_LD', 'age_class', 'temp_mean_c', 'year_c') %in% names(d)))
d$age_class <- factor(d$age_class, levels = c('YOUNG', 'OLD'))
X <- model.matrix(~ rel_LD + age_class + temp_mean_c + year_c, d)
stopifnot(!anyNA(X), qr(X)$rank == ncol(X))
Z <- X[, colnames(X) != '(Intercept)', drop = FALSE]
vif <- vapply(seq_len(ncol(Z)), function(j) {
  fit <- lm(Z[, j] ~ Z[, -j, drop = FALSE])
  1 / (1 - summary(fit)$r.squared)
}, numeric(1))
annual <- unique(d[, c('year_f', 'temp_mean_c', 'year_c')])
stopifnot(nrow(annual) == 39L)
after <- unname(tools::md5sum(model_path))
stopifnot(identical(before, after))
cat('FIXED_DESIGN_VIF_CSV\n')
write.csv(data.frame(predictor = colnames(Z), VIF = vif), stdout(), row.names = FALSE)
cat('METADATA_CSV\n')
write.csv(data.frame(n = nrow(d), females = length(unique(d$female_id)),
  seasons = nrow(annual), min_year = min(as.integer(as.character(d$year_f))),
  max_year = max(as.integer(as.character(d$year_f))), design_rank = qr(X)$rank,
  annual_r_temperature_year = cor(annual$temp_mean_c, annual$year_c),
  brood_r_temperature_year = cor(d$temp_mean_c, d$year_c),
  brood_r_date_temperature = cor(d$rel_LD, d$temp_mean_c),
  brood_r_date_year = cor(d$rel_LD, d$year_c),
  max_VIF = max(vif), model_md5_unchanged = identical(before, after)),
  stdout(), row.names = FALSE)
cat('MODEL_FORMULA\n')
print(m$formula)
cat('SESSION\n')
print(sessionInfo())
# Both location and log precision have this additive fixed-effect design.
# The constant extra-zero component has no slope VIF. These are unweighted
# design-matrix VIFs, not posterior-covariance or mixed-model precision measures.
