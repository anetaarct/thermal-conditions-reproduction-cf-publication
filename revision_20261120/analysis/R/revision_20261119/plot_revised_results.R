root <- '.'
out <- 'figures'
dir.create(out, recursive = TRUE, showWarnings = FALSE)
baseline <- read.csv(file.path(root, 'tables/reproduced/baseline_predictions/curves.csv'))
within <- read.csv(file.path(root, 'tables/within_between_20261026_clean_20261119/prediction_curves.csv'))

draw <- function(data, specs, name, ncol = 2) {
  panel_rows <- ceiling(length(specs) / ncol)
  pdf(file.path(out, name), width = 9, height = 3.5 * panel_rows, useDingbats = FALSE)
  par(mfrow = c(panel_rows, ncol), mar = c(4.1, 4.8, 2, 1), mgp = c(2.8, .8, 0), las = 1)
  for (i in seq_along(specs)) {
    spec <- specs[[i]]
    z <- subset(data, trait == spec$trait & predictor == spec$predictor & component == spec$component)
    stopifnot(nrow(z) > 0)
    percentage <- spec$trait == 'fledging' && spec$component == 'mean'
    mult <- if (percentage) 100 else 1
    discrete <- spec$predictor == 'age_class'
    x <- if (discrete) match(z$x, c('YOUNG', 'OLD')) else as.numeric(z$x)
    stopifnot(!anyNA(x), all(is.finite(z$estimate)), all(is.finite(z$lower)), all(is.finite(z$upper)))
    ord <- order(x); x <- x[ord]; z <- z[ord, ]
    ylabel <- if (spec$component == 'mean') {
      if (spec$trait == 'clutch') 'Expected clutch size (eggs)' else 'Expected egg-to-fledging success (%)'
    } else if (spec$trait == 'clutch') 'Conditional clutch-size SD (eggs)' else 'Beta-binomial precision'
    xlabel <- switch(spec$predictor, temp_mean_c = 'Spring temperature deviation (degrees C)',
      temp_within = 'Within-female temperature deviation (degrees C)',
      temp_between = 'Female mean spring temperature (degrees C)',
      rel_LD = 'Relative laying date (days)', age_class = 'Female age class', year = 'Calendar year')
    plot(x, z$estimate * mult, type = 'n', ylim = range(z$lower, z$upper) * mult,
      xlim = if (discrete) c(.5, 2.5) else range(x), xlab = xlabel, ylab = ylabel,
      bty = 'l', xaxt = if (discrete) 'n' else 's')
    if (discrete) {
      axis(1, at = 1:2, labels = c('Yearling', 'Older'))
      arrows(x, z$lower * mult, x, z$upper * mult, angle = 90, code = 3, length = .05, col = '#0072B2')
      points(x, z$estimate * mult, pch = 16, col = '#0072B2')
    } else {
      polygon(c(x, rev(x)), c(z$lower, rev(z$upper)) * mult, col = adjustcolor('#0072B2', .18), border = NA)
      lines(x, z$estimate * mult, col = '#0072B2', lwd = 2)
    }
    title(LETTERS[i], adj = 0)
  }
  dev.off()
}

spec <- function(trait, predictor, component) list(trait = trait, predictor = predictor, component = component)
four <- function(predictor) list(spec('clutch', predictor, 'mean'), spec('clutch', predictor, 'consistency'),
  spec('fledging', predictor, 'mean'), spec('fledging', predictor, 'consistency'))
draw(baseline, list(spec('clutch', 'temp_mean_c', 'mean'), spec('fledging', 'temp_mean_c', 'mean')), 'Fig_temperature_main_20261119.pdf')
draw(baseline, four('temp_mean_c'), 'Fig_temperature_full_20261119.pdf')
draw(baseline, four('rel_LD'), 'Fig_laying_date_20261119.pdf')
draw(baseline, four('age_class'), 'Fig_age_class_20261119.pdf')
for (trait in c('clutch', 'fledging')) {
  panels <- list(spec(trait, 'temp_within', 'mean'), spec(trait, 'temp_between', 'mean'),
    spec(trait, 'temp_within', 'consistency'), spec(trait, 'temp_between', 'consistency'))
  draw(within, panels, paste0('Fig_within_between_', trait, '_20261119.pdf'))
}
temporal <- list()
for (trait in c('clutch', 'fledging')) {
  components <- if (trait == 'clutch') c(mean = 'mean', sd = 'consistency') else c(success = 'mean', phi = 'consistency')
  for (component in names(components)) {
    z <- read.csv(file.path(root, 'tables/temporal_clean_20261119', paste0(trait, '_year_', component, '.csv')))
    temporal[[paste(trait, component)]] <- data.frame(trait, predictor = 'year', component = components[[component]],
      x = z$x, estimate = z$estimate, lower = z$CI_low, upper = z$CI_high)
  }
}
draw(do.call(rbind, temporal), four('year'), 'Fig_temporal_without_temperature_20261119.pdf')
cat('All revised reproductive figures created\n')
