# Validate the public model-ready records against final published sample summaries.
checksums<-read.csv('data_derived/analysis/checksums.csv',stringsAsFactors=FALSE)
actual<-tools::md5sum(file.path('data_derived/analysis',checksums$filename))
stopifnot(!anyNA(actual),identical(unname(actual),checksums$md5))
metadata <- read.csv('tables/within_between_20261026/sample_metadata.csv')
reference <- read.csv('tables/temporal_models_20261022/observed_annual.csv')
checks <- list()
for(trait in c('clutch','fledging')) {
  d <- read.csv(paste0('data_derived/analysis/',trait,'_model_data.csv'))
  m <- metadata[metadata$trait==trait,]
  stopifnot(nrow(d)==m$records,length(unique(d$female_id))==m$females,
            length(unique(d$year))==m$seasons,!anyNA(d),
            all(d$clutch %in% 2:9),all(d$age_class %in% c('YOUNG','OLD')),
            all(as.integer(d$year_f)==d$year))
  if(trait=='fledging') stopifnot(all(d$fledged>=0 & d$fledged<=d$clutch))
  observed <- do.call(rbind,lapply(split(d,d$year),function(x) {
    y <- if(trait=='clutch') x$clutch else x$fledged/x$clutch
    data.frame(year=x$year[1],n=nrow(x),observed_mean=mean(y),observed_SD=sd(y),
      older_fraction=mean(x$age_class=='OLD'),mean_relative_date=mean(x$rel_LD),
      mean_temperature=mean(x$temp_mean_c))
  }))
  ref <- reference[reference$trait==trait,]
  observed <- observed[match(ref$year,observed$year),]
  for(column in setdiff(names(observed),'year')) {
    delta <- max(abs(observed[[column]]-ref[[column]]))
    stopifnot(delta<1e-9)
    checks[[paste(trait,column)]] <- data.frame(trait,check=column,max_absolute_difference=delta)
  }
  groups <- split(d$temp_mean,d$female_id)
  varied <- vapply(groups,function(x) length(unique(round(x,10)))>1,logical(1))
  stopifnot(sum(varied)==m$females_varied_temperature,
            sum(lengths(groups)[varied])==m$records_in_varied_temperature_females,
            max(abs(d$temp_mean-d$temp_mean_c-m$temperature_center))<1e-8,
            max(abs(d$year-d$year_c-m$year_center))<1e-8)
}
dir.create('tables/reproduced',recursive=TRUE,showWarnings=FALSE)
write.csv(do.call(rbind,checks),'tables/reproduced/input_validation.csv',row.names=FALSE)
print(do.call(rbind,checks),row.names=FALSE)
cat('Public model data match final sample summaries.\n')
