.libPaths(c(path.expand('~/R/cf_model_comparison_20261013_R46'),path.expand('~/R/x86_64-pc-linux-gnu-library/4.5'),.libPaths()))
library(glmmTMB)
rows<-list()
for(kind in c('baseline','within_between')){
 dir<-file.path('remaining_robustness_20261120',paste('scale_year_clutch',kind,sep='_'))
 reference<-readRDS(paste0('models/clean_20261119/clutch_',kind,'_selected.rds'))
 for(method in c('NLMINB','BFGS')){
  file<-file.path(dir,paste0('scale_year_clutch_',kind,'_',method,'.rds'))
  if(!file.exists(file))next
  m<-readRDS(file)
  cf<-fixef(m);orig<-fixef(reference)
  rows[[paste(kind,method)]]<-data.frame(model=kind,optimizer=method,convergence=m$fit$convergence,pdHess=isTRUE(m$sdr$pdHess),year_scale_SD=unname(attr(VarCorr(m)$disp$year_f,'stddev')),logLik=as.numeric(logLik(m)),objective=m$fit$objective,reference_logLik=as.numeric(logLik(reference)),logLik_gain=-m$fit$objective-as.numeric(logLik(reference)),delta_AIC=AIC(m)-AIC(reference),max_abs_cond_change=max(abs(cf$cond-orig$cond)),max_abs_disp_change=max(abs(cf$disp-orig$disp)))
 }
}
write.csv(do.call(rbind,rows),'remaining_robustness_20261120/scale_boundary_audit.csv',row.names=FALSE);print(do.call(rbind,rows))
