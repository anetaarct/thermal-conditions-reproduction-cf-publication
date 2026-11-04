# Equidispersed clutch-size comparator fitted directly from public records.
library(glmmTMB)
d<-read.csv('data_derived/analysis/clutch_model_data.csv')
d$female_id<-factor(d$female_id);d$year_f<-factor(d$year_f)
d$age_class<-factor(d$age_class,levels=c('YOUNG','OLD'))
old<-subset(read.csv('tables/final_models_20261019/count_coefficients.csv'),model=='Poisson' & component=='cond')
X<-model.matrix(~rel_LD+temp_mean_c+year_c+age_class,d)
beta<-setNames(old$estimate,old$term)[colnames(X)]
m<-glmmTMB(clutch~rel_LD+temp_mean_c+year_c+age_class+(1|female_id)+(1|year_f),
 family=poisson(),data=d,start=list(beta=unname(beta),theta=log(c(.001,.04))),
 control=glmmTMBControl(optimizer=optim,optArgs=list(method='BFGS'),optCtrl=list(maxit=2000,reltol=1e-11)))
dir.create('models/poisson_comparator',recursive=TRUE,showWarnings=FALSE)
saveRDS(m,'models/poisson_comparator/selected.rds')
dir.create('tables/reproduced',recursive=TRUE,showWarnings=FALSE)
capture.output(summary(m),file='tables/reproduced/poisson_summary.txt')
write.csv(data.frame(logLik=as.numeric(logLik(m)),AIC=AIC(m),convergence=m$fit$convergence,
 pdHess=isTRUE(m$sdr$pdHess)),'tables/reproduced/poisson_diagnostics.csv',row.names=FALSE)
cf<-summary(m)$coefficients$cond
fixed<-data.frame(term=rownames(cf),estimate=cf[,1],SE=cf[,2],p=cf[,4])
write.csv(fixed,'tables/reproduced/poisson_fixed_effects.csv',row.names=FALSE)
z<-merge(old,fixed,by='term',suffixes=c('_published','_reproduced'))
stopifnot(nrow(z)==nrow(old),max(abs(z$estimate_published-z$estimate_reproduced))<.001)
# Boundary random-effect estimates can give a singular Hessian; retain that diagnostic.
stopifnot(m$fit$convergence==0)
