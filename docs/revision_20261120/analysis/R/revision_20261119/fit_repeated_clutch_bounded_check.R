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
 z<-optim(par,fn,gr,method='L-BFGS-B',lower=lower,upper=upper,control=list(maxit=1500,factr=1e5,pgtol=1e-7,parscale=scale,trace=1,REPORT=1))
 list(par=z$par,objective=z$value,convergence=z$convergence,message=z$message)
}
a<-glmmTMB(formula(reference),dispformula=formula(reference,component='disp'),family=compois(),data=d,start=st,control=glmmTMBControl(optimizer=optimizer,parallel=1))
saveRDS(a,file.path(out,'bounded_start.rds'))
cat('BOUNDED COMPLETE',as.numeric(logLik(a)),a$fit$convergence,a$sdr$pdHess,'boundary distance',min(pmin(a$fit$par-lower,upper-a$fit$par)/(upper-lower)),'\n');flush.console()
capture.output(summary(a),file=file.path(out,'bounded_summary.txt'))
V<-a$sdr$cov.fixed;e<-eigen(V,symmetric=TRUE);stopifnot(min(e$values)>0);S<-e$vectors%*%diag(sqrt(e$values))
p<-a$obj$env$parList(par=a$fit$parfull);st<-p[intersect(c('beta','betadisp','theta','b'),names(p))]
independent<-function(par,fn,gr,control=list(),...){
 inner<-environment(fn)$inner.control;inner$maxit<-1000L;inner$grad.tol<-1e-9;environment(fn)$inner.control<-inner
 ff<-function(q)fn(drop(par+S%*%q));gg<-function(q)drop(crossprod(S,gr(drop(par+S%*%q))))
 z<-nlminb(rep(0,length(par)),ff,gg,control=list(iter.max=1500,eval.max=3000,rel.tol=1e-10,trace=1))
 list(par=drop(par+S%*%z$par),objective=z$objective,convergence=z$convergence,message=z$message)
}
b<-update(a,data=d,start=st,control=glmmTMBControl(optimizer=independent,parallel=1))
saveRDS(b,file.path(out,'unbounded_independent.rds'))
dg<-do.call(rbind,lapply(list(LBFGSB=a,NLMINB=b),function(x)data.frame(convergence=x$fit$convergence,pdHess=x$sdr$pdHess,logLik=as.numeric(logLik(x)),AIC=AIC(x),scaled_gradient=max(abs(x$sdr$gradient.fixed*sqrt(diag(x$sdr$cov.fixed)))))))
write.csv(dg,file.path(out,'bounded_independent_diagnostics.csv'));print(dg)
stopifnot(all(dg$convergence==0),all(dg$pdHess),all(dg$scaled_gradient<.005),diff(range(dg$logLik))<.001,min(pmin(a$fit$par-lower,upper-a$fit$par)/(upper-lower))>.001)
saveRDS(b,file.path(out,'repeated_females_clutch_within_between_selected.rds'))
code<-readLines('export_sensitivity_contrasts.R');eval(parse(text=code[seq_len(grep('^age<-list',code)[1]-1)]))
effects<-do.call(rbind,lapply(c('temp_within','temp_between'),function(pr)cbind(data.frame(check='repeated_females',trait='clutch',model='within_between',predictor=pr),evaluate_contrast(b,'clutch',pr,20261120)[,c('mean_change','lower','upper','consistency_ratio','ratio_lower','ratio_upper')])))
write.csv(effects,file.path(out,'effects.csv'),row.names=FALSE);print(effects)
cat('BOUNDED AND UNBOUNDED CHECK COMPLETE\n')
