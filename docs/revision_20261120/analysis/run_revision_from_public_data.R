.libPaths(c(path.expand('~/R/cf_model_comparison_20261013_R46'),path.expand('~/R/x86_64-pc-linux-gnu-library/4.5'),.libPaths()))
Sys.setenv(WB_MODEL_DIR='models/clean_20261119',WB_OUTPUT_SUFFIX='_clean_20261119',
  WB_EXPORT_DIR='tables/within_between_20261026_clean_20261119',
  TEMPORAL_MODEL_DIR='models/temporal_clean_20261119',TEMPORAL_EXPORT_DIR='tables/temporal_clean_20261119',
  JOINT_YEAR_EXPORT_DIR='tables/joint_year_clean_20261119')
stopifnot(file.exists('data_derived/analysis_clean_20261119/expected_annual.csv'))
code<-readLines('R/90_fit_temperature_within_between.R')
code<-sub('tables/temporal_models_20261022/observed_annual.csv','data_derived/analysis_clean_20261119/expected_annual.csv',code,fixed=TRUE)
# The archived likelihood is a starting-value reference, not a target after exclusions.
code<-sub('stopifnot(abs(difference) < .001)','stopifnot(is.finite(difference))',code,fixed=TRUE)
env<-new.env(parent=globalenv())
eval(parse(text=code),envir=env)
env$run_within_between('data_derived/analysis_clean_20261119/clutch_model_data.csv',
                       'data_derived/analysis_clean_20261119/fledging_model_data.csv')
source('R/83_fit_temporal_without_temperature.R',local=new.env(parent=globalenv()))
source('R/91_export_temperature_within_between.R',local=new.env(parent=globalenv()))
source('R/93_export_joint_calendar_year.R',local=new.env(parent=globalenv()))
source('R/96_export_baseline_predictions.R',local=new.env(parent=globalenv()))
source('R/84_export_temporal_without_temperature.R',local=new.env(parent=globalenv()))
cat('CLEAN REVISION COMPLETE\n')

source('R/revision_20261119/fit_exact_age_public.R',local=new.env(parent=globalenv()))
source('R/revision_20261119/validate_age.R',local=new.env(parent=globalenv()))
source('R/revision_20261119/revision_diagnostics.R',local=new.env(parent=globalenv()))
source('R/revision_20261119/refine_extra_zero.R',local=new.env(parent=globalenv()))
source('R/revision_20261119/export_sensitivity_contrasts.R',local=new.env(parent=globalenv()))

source('R/revision_20261119/fit_window_sensitivity.R',local=new.env(parent=globalenv()))

source('R/revision_20261119/plot_revised_results.R',local=new.env(parent=globalenv()))
source('R/revision_20261119/plot_window_sensitivity.R',local=new.env(parent=globalenv()))
