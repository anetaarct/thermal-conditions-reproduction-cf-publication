.libPaths(c(path.expand('~/R/cf_model_comparison_20261013_R46'),path.expand('~/R/x86_64-pc-linux-gnu-library/4.5'),.libPaths()))
library(glmmTMB)
library(MASS)
evaluate_contrast <- function(m,trait,predictor,seed) {
 cf<-fixef(m);nc<-length(cf$cond);nd<-length(cf$disp)
 names_all<-c(names(cf$cond),paste0('disp~',names(cf$disp)),if(length(cf$zi))paste0('zi~',names(cf$zi)))
 V<-vcov(m,full=TRUE)[names_all,names_all,drop=FALSE]
 set.seed(seed)
 draws<-mvrnorm(6000,c(unname(cf$cond),unname(cf$disp),unname(cf$zi)),V)
 make_x<-function(b) {
   x<-setNames(rep(0,length(b)),names(b));x['(Intercept)']<-1
   frame<-model.frame(m)
   if('age_c' %in% names(b))x['age_c']<-min(frame$age_c)
   if('age_c2' %in% names(b))x['age_c2']<-x['age_c']^2
   x
 }
 X<-rbind(make_x(cf$cond),make_x(cf$cond));X[2,predictor]<-1
 Z<-rbind(make_x(cf$disp),make_x(cf$disp));Z[2,predictor]<-1
 eta<-drop(X%*%cf$cond);etas<-X%*%t(draws[,seq_len(nc),drop=FALSE])
 phi<-exp(drop(Z%*%cf$disp));phis<-exp(Z%*%t(draws[,nc+seq_len(nd),drop=FALSE]))
 if(trait=='clutch') {
   mu<-exp(eta);mus<-exp(etas)
   sd<-sqrt(compois()$variance(mu,phi))
   sds<-matrix(sqrt(compois()$variance(as.vector(mus),as.vector(phis))),2,nrow(draws))
   ratio<-sds[2,]/sds[1,];sd_delta<-sds[2,]-sds[1,]
   data.frame(mean_first=mu[1],mean_last=mu[2],mean_change=diff(mu),
     lower=quantile(mus[2,]-mus[1,],.025),upper=quantile(mus[2,]-mus[1,],.975),
     SD_first=sd[1],SD_last=sd[2],SD_change=diff(sd),SD_change_lower=quantile(sd_delta,.025),SD_change_upper=quantile(sd_delta,.975),
     consistency_ratio=sd[2]/sd[1],ratio_lower=quantile(ratio,.025),ratio_upper=quantile(ratio,.975))
 } else {
   W<-rbind(make_x(cf$zi),make_x(cf$zi))
   if(predictor %in% colnames(W))W[2,predictor]<-1
   pi<-plogis(drop(W%*%cf$zi));pis<-plogis(W%*%t(draws[,nc+nd+seq_along(cf$zi),drop=FALSE]))
   mu<-plogis(eta);mus<-plogis(etas)
   success<-(1-pi)*mu;successes<-(1-pis)*mus
   zero<-pi+(1-pi)*exp(lbeta(mu*phi,(1-mu)*phi+6)-lbeta(mu*phi,(1-mu)*phi))
   zeros<-pis+(1-pis)*exp(lbeta(mus*phis,(1-mus)*phis+6)-lbeta(mus*phis,(1-mus)*phis))
   data.frame(mean_first=success[1],mean_last=success[2],mean_change=diff(success),
     lower=quantile(successes[2,]-successes[1,],.025),upper=quantile(successes[2,]-successes[1,],.975),
     failure_first=zero[1],failure_last=zero[2],failure_change=diff(zero),
     failure_change_lower=quantile(zeros[2,]-zeros[1,],.025),failure_change_upper=quantile(zeros[2,]-zeros[1,],.975),
     consistency_ratio=phi[2]/phi[1],ratio_lower=quantile(phis[2,]/phis[1,],.025),ratio_upper=quantile(phis[2,]/phis[1,],.975))
 }
}
age<-list()
for(trait in c('clutch','fledging')) for(kind in c('age_class','exact_age_linear','exact_age_quadratic')) {
 path<-file.path('age_sensitivity_20261119',paste0(trait,'_',kind,'_refined.rds'))
 if(!file.exists(path))next
 m<-readRDS(path)
 if(!isTRUE(m$sdr$pdHess) || m$fit$convergence!=0)next
 z<-evaluate_contrast(m,trait,'temp_within',20261119)
 # Keep common fields across outcome distributions.
 age[[paste(trait,kind)]]<-cbind(data.frame(trait,model=kind),z[,c('mean_change','lower','upper','consistency_ratio','ratio_lower','ratio_upper')])
}
write.csv(do.call(rbind,age),'age_sensitivity_20261119/reference_contrasts.csv',row.names=FALSE)
extra<-list()
for(kind in c('year_only','temperature')) {
 path<-if(kind=='year_only')'models/clean_20261119/fledging_baseline_selected.rds' else 'tables/curation_checks_20261119/extra_zero_temperature_refined.rds'
 if(!file.exists(path))next
 m<-readRDS(path)
 if(!isTRUE(m$sdr$pdHess) || m$fit$convergence!=0)next
 extra[[kind]]<-cbind(data.frame(model=kind),evaluate_contrast(m,'fledging','temp_mean_c',20261120))
}
write.csv(do.call(rbind,extra),'tables/curation_checks_20261119/full_mixture_contrasts.csv',row.names=FALSE)
clutch<-readRDS('models/clean_20261119/clutch_within_between_selected.rds')
sd<-do.call(rbind,lapply(c('temp_within','temp_between'),function(p)cbind(data.frame(predictor=p),evaluate_contrast(clutch,'clutch',p,20261121))))
write.csv(sd,'tables/curation_checks_20261119/clutch_SD_contrasts.csv',row.names=FALSE)
random<-stats<-list()
for(trait in c('clutch','fledging')) {
 m<-readRDS(file.path('models/clean_20261119',paste0(trait,'_baseline_selected.rds')))
 d<-readRDS(file.path('models/clean_20261119',paste0(trait,'_analysis_data.rds')))
 re<-VarCorr(m)
 random[[trait]]<-do.call(rbind,lapply(names(re),function(component)do.call(rbind,lapply(names(re[[component]]),function(group)
   data.frame(trait,component,group,SD=unname(attr(re[[component]][[group]],'stddev')))))))
 stats[[trait]]<-data.frame(trait,n=nrow(d),females=length(unique(d$female_id)),seasons=length(unique(d$year_f)),
    mean_clutch=mean(d$clutch),variance_clutch=var(d$clutch),
    complete_failures=if(trait=='fledging')sum(d$fledged==0) else NA,
    mean_success=if(trait=='fledging')mean(d$fledged/d$clutch) else NA)
}
write.csv(do.call(rbind,random),'tables/curation_checks_20261119/baseline_random_SD.csv',row.names=FALSE)
write.csv(do.call(rbind,stats),'tables/curation_checks_20261119/sample_summary.csv',row.names=FALSE)
for(trait in c('clutch','fledging')) {
 p<-file.path('age_sensitivity_20261119',paste0(trait,'_exact_age_linear_refined.rds'))
 if(!file.exists(p))next
 m<-readRDS(p);f<-model.frame(m)
 age<-as.integer(round(f$age_c-min(f$age_c)+1))
 d<-read.csv(file.path('data_derived/analysis_clean_20261119',paste0(trait,'_model_data.csv')))
 keys<-paste(d$female_id,d$year,sep='|');modelkeys<-paste(f$female_id,f$year_f,sep='|')
 stopifnot(!anyDuplicated(keys),all(modelkeys %in% keys))
 matched<-d[match(modelkeys,keys),];matched$exact_age<-age
 write.csv(matched,file.path('age_sensitivity_20261119',paste0(trait,'_known_age_model_data.csv')),row.names=FALSE)
}
cat('SENSITIVITY CONTRASTS COMPLETE\n')
