.libPaths(c(path.expand('~/R/cf_model_comparison_20261013_R46'),path.expand('~/R/x86_64-pc-linux-gnu-library/4.5'),.libPaths()))
library(glmmTMB)
out<-'tables/curation_checks_20261119'
d<-readRDS('models/clean_20261119/fledging_analysis_data.rds')
baseline<-readRDS('models/clean_20261119/fledging_baseline_selected.rds')
pl<-baseline$obj$env$parList(par=baseline$fit$parfull)
start<-pl[intersect(c('beta','betazi','betadisp','theta','thetadisp'),names(pl))]
start<-start[lengths(start)>0]
m<-glmmTMB(cbind(fledged,clutch-fledged)~rel_LD+age_class+temp_mean_c+year_c+(1|female_id)+(1|year_f),
 dispformula=~rel_LD+age_class+temp_mean_c+year_c+(1|year_f),ziformula=~1,
 family=betabinomial(),data=d,start=start,
 control=glmmTMBControl(optimizer=optim,optArgs=list(method='BFGS'),optCtrl=list(maxit=5000,reltol=1e-12),parallel=4))
saveRDS(m,file.path(out,'extra_zero_constant_refined.rds'))
a<-readRDS(file.path(out,'extra_zero_temperature.rds'))
V<-a$sdr$cov.fixed
optimizer<-function(par,fn,gr,control=list(),...) {
 ee<-eigen((V+t(V))/2,symmetric=TRUE);stopifnot(min(ee$values)>0)
 S<-ee$vectors%*%diag(sqrt(ee$values))
 z<-optim(rep(0,length(par)),function(q)fn(drop(par+S%*%q)),function(q)drop(crossprod(S,gr(drop(par+S%*%q)))),method='BFGS',control=list(maxit=2000,reltol=1e-12))
 list(par=drop(par+S%*%z$par),objective=z$value,convergence=z$convergence,message=z$message)
}
pl<-a$obj$env$parList(par=a$fit$parfull)
st<-pl[intersect(c('beta','betazi','betadisp','theta','thetadisp','thetazi'),names(pl))]
st<-st[lengths(st)>0]
b<-glmmTMB(formula(a),dispformula=formula(a,component='disp'),ziformula=formula(a,component='zi'),
 family=betabinomial(),data=d,start=st,control=glmmTMBControl(optimizer=optimizer,parallel=4))
saveRDS(b,file.path(out,'extra_zero_temperature_refined.rds'))
models<-list(year_only=baseline,constant=m,temperature=b)
diag<-do.call(rbind,lapply(names(models),function(kind){x<-models[[kind]];data.frame(model=kind,n=nobs(x),
 convergence=x$fit$convergence,pdHess=isTRUE(x$sdr$pdHess),logLik=as.numeric(logLik(x)),AIC=AIC(x),
 scaled_gradient=max(abs(x$sdr$gradient.fixed*sqrt(diag(x$sdr$cov.fixed)))))}))
write.csv(diag,file.path(out,'extra_zero_validated_comparison.csv'),row.names=FALSE)
co<-do.call(rbind,lapply(names(models),function(kind){x<-summary(models[[kind]])$coefficients;do.call(rbind,lapply(names(x),function(component){v<-x[[component]];if(!length(v))return(NULL)
 data.frame(model=kind,component,term=rownames(v),estimate=v[,1],SE=v[,2],lower=v[,1]-qnorm(.975)*v[,2],upper=v[,1]+qnorm(.975)*v[,2])}))}))
write.csv(co,file.path(out,'extra_zero_validated_fixed.csv'),row.names=FALSE)
print(diag);cat('Temperature-model logLik refinement difference',as.numeric(logLik(b)-logLik(a)),'\n')
cat('EXTRA ZERO VALIDATION COMPLETE\n')
