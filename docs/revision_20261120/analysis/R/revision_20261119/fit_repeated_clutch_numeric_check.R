.libPaths(c(path.expand('~/R/cf_model_comparison_20261013_R46'),path.expand('~/R/x86_64-pc-linux-gnu-library/4.5'),.libPaths()))
library(glmmTMB)
out<-'remaining_robustness_20261120/repeated_females_clutch_within_between'
m<-readRDS('remaining_robustness_20261120/repeated_females_clutch_baseline/repeated_females_clutch_baseline_selected.rds')
d<-readRDS('models/clean_20261119/clutch_analysis_data.rds');counts<-table(d$female_id);d<-droplevels(d[as.character(d$female_id)%in%names(counts[counts>1]),])
reference<-readRDS('models/clean_20261119/clutch_within_between_selected.rds')
cols<-colnames(model.matrix(~rel_LD+temp_within+temp_between+year_c+age_class,d))
p<-m$obj$env$parList(par=m$fit$parfull);st<-p[intersect(c('beta','betadisp','theta','b'),names(p))]
for(component in c('cond','disp')){b<-fixef(m)[[component]];b['temp_within']<-b['temp_mean_c'];b['temp_between']<-b['temp_mean_c'];st[[if(component=='cond')'beta' else 'betadisp']]<-unname(b[cols])}
lower<-upper<-scale<-NULL
optimizer<-function(par,fn,gr,control=list(),...){
 inner<-environment(fn)$inner.control;inner$maxit<-1000L;inner$grad.tol<-1e-9;environment(fn)$inner.control<-inner
 lower<<-rep(-10,length(par));upper<<-rep(0,length(par));scale<<-rep(1,length(par))
 for(g in c('beta','betadisp')){
  ix<-which(names(par)==g);lower[ix]<<-if(g=='beta')-.5 else -1;upper[ix]<<-if(g=='beta').5 else 1
  lower[ix[cols=='(Intercept)']]<<-if(g=='beta')1.5 else -6
  upper[ix[cols=='(Intercept)']]<<-if(g=='beta')2.1 else -1
  lower[ix[cols=='year_c']]<<-if(g=='beta')-.02 else -.1
  upper[ix[cols=='year_c']]<<-if(g=='beta').02 else .1
  scale[ix]<<-if(g=='beta')c(.004,.0003,.004,.004,.0004,.003) else c(.04,.003,.04,.03,.002,.04)
 }
 z<-optim(par,fn,NULL,method='L-BFGS-B',lower=lower,upper=upper,control=list(maxit=100,factr=1e6,pgtol=1e-7,parscale=scale,trace=1,REPORT=1))
 list(par=z$par,objective=z$value,convergence=z$convergence,message=z$message)
}
a<-glmmTMB(formula(reference),dispformula=formula(reference,component='disp'),family=compois(),data=d,start=st,control=glmmTMBControl(optimizer=optimizer,parallel=1))
saveRDS(a,file.path(out,'numeric_gradient_start.rds'))
cat('BOUNDED COMPLETE',as.numeric(logLik(a)),a$fit$convergence,a$sdr$pdHess,'boundary distance',min(pmin(a$fit$par-lower,upper-a$fit$par)/(upper-lower)),'\n');flush.console()
capture.output(summary(a),file=file.path(out,'numeric_gradient_summary.txt'))

write.csv(data.frame(convergence=a$fit$convergence,pdHess=a$sdr$pdHess,logLik=as.numeric(logLik(a)),AIC=AIC(a),scaled_gradient=max(abs(a$sdr$gradient.fixed*sqrt(diag(a$sdr$cov.fixed))))),file.path(out,'numeric_gradient_diagnostics.csv'),row.names=FALSE)
