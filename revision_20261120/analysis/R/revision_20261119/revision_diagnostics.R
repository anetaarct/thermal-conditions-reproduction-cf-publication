.libPaths(c(path.expand('~/R/cf_model_comparison_20261013_R46'),path.expand('~/R/x86_64-pc-linux-gnu-library/4.5'),.libPaths()))
Sys.setenv(WB_MODEL_DIR='models/clean_20261119')
library(glmmTMB)
out <- 'tables/curation_checks_20261119'
dir.create(out,recursive=TRUE,showWarnings=FALSE)
for(trait in c('clutch','fledging')) {
  f<-file.path('models/clean_20261119',paste0(trait,'_baseline_selected.rds'))
  while(!file.exists(f))Sys.sleep(15)
}
source('R/98_baseline_simulation_checks.R',local=new.env(parent=globalenv()))
checks<-vifs<-list()
for(trait in c('clutch','fledging')) {
  m<-readRDS(file.path('models/clean_20261119',paste0(trait,'_baseline_selected.rds')))
  d<-readRDS(file.path('models/clean_20261119',paste0(trait,'_analysis_data.rds')))
  X<-model.matrix(~rel_LD+temp_mean_c+year_c+age_class,d)[,-1,drop=FALSE]
  v<-vapply(seq_len(ncol(X)),function(j)1/(1-summary(lm(X[,j]~X[,-j,drop=FALSE]))$r.squared),numeric(1))
  vifs[[trait]]<-data.frame(trait,term=colnames(X),VIF=v)
  set.seed(20261119+match(trait,c('clutch','fledging')))
  simulated<-simulate(m,nsim=300)
  draws<-vapply(simulated,function(x)if(is.matrix(x))x[,1] else as.numeric(x),numeric(nrow(d)))
  observed<-if(trait=='clutch')d$clutch else d$fledged/d$clutch
  if(trait=='fledging')draws<-draws/d$clutch
  stats<-list(mean=colMeans(draws),variance=apply(draws,2,var))
  obs<-c(mean=mean(observed),variance=var(observed))
  for(s in names(stats)) checks[[paste(trait,s)]]<-data.frame(trait,target='new_group_plugin',statistic=s,observed=obs[s],
    predicted=mean(stats[[s]]),lower=quantile(stats[[s]],.025),upper=quantile(stats[[s]],.975))
  write.csv(do.call(rbind,checks),file.path(out,'new_group_simulations.csv'),row.names=FALSE)
  write.csv(do.call(rbind,vifs),file.path(out,'baseline_vif.csv'),row.names=FALSE)
}
# Matched, corrected-sample Poisson comparator.
d<-readRDS('models/clean_20261119/clutch_analysis_data.rds')
poisson_file<-file.path(out,'poisson.rds')
if(file.exists(poisson_file))p<-readRDS(poisson_file) else {
 p<-glmmTMB(clutch~rel_LD+temp_mean_c+year_c+age_class+(1|female_id)+(1|year_f),
   family=poisson(),data=d,control=glmmTMBControl(optCtrl=list(iter.max=3000,eval.max=6000),parallel=4))
 saveRDS(p,poisson_file)
}
write.csv(data.frame(logLik=as.numeric(logLik(p)),AIC=AIC(p),convergence=p$fit$convergence,pdHess=isTRUE(p$sdr$pdHess)),
          file.path(out,'poisson_diagnostics.csv'),row.names=FALSE)
co<-summary(p)$coefficients$cond
write.csv(data.frame(term=rownames(co),estimate=co[,1],SE=co[,2],lower=co[,1]-qnorm(.975)*co[,2],upper=co[,1]+qnorm(.975)*co[,2]),file.path(out,'poisson_fixed.csv'),row.names=FALSE)

# Separate environmental-slope check for the additional-zero component.
d<-readRDS('models/clean_20261119/fledging_analysis_data.rds')
models<-list(year_only=readRDS('models/clean_20261119/fledging_baseline_selected.rds'))
for(kind in c('constant','temperature')) {
 file<-file.path(out,paste0('extra_zero_',kind,'.rds'))
 if(file.exists(file))m<-readRDS(file) else {
  m<-glmmTMB(cbind(fledged,clutch-fledged)~rel_LD+age_class+temp_mean_c+year_c+(1|female_id)+(1|year_f),
    dispformula=~rel_LD+age_class+temp_mean_c+year_c+(1|year_f),
    ziformula=if(kind=='constant')~1 else ~1+temp_mean_c+(1|year_f),family=betabinomial(),data=d,
    control=glmmTMBControl(optCtrl=list(iter.max=3000,eval.max=6000),parallel=4))
  saveRDS(m,file)
 }
 models[[kind]]<-m
}
diag<-do.call(rbind,lapply(names(models),function(kind){m<-models[[kind]];data.frame(model=kind,n=nobs(m),
  convergence=m$fit$convergence,pdHess=isTRUE(m$sdr$pdHess),logLik=as.numeric(logLik(m)),AIC=AIC(m),
  scaled_gradient=max(abs(m$sdr$gradient.fixed*sqrt(diag(m$sdr$cov.fixed)))))}))
write.csv(diag,file.path(out,'extra_zero_model_comparison.csv'),row.names=FALSE)
co<-do.call(rbind,lapply(names(models),function(kind){m<-models[[kind]];do.call(rbind,lapply(names(fixef(m)),function(component){
 x<-summary(m)$coefficients[[component]];if(!length(x))return(NULL)
 data.frame(model=kind,component,term=rownames(x),estimate=x[,1],SE=x[,2],lower=x[,1]-qnorm(.975)*x[,2],upper=x[,1]+qnorm(.975)*x[,2])}))}))
write.csv(co,file.path(out,'extra_zero_fixed.csv'),row.names=FALSE)
capture.output(sessionInfo(),file=file.path(out,'session_info.txt'))
cat('ADDITIONAL DIAGNOSTICS COMPLETE\n')
