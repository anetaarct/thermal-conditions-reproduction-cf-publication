.libPaths(c(path.expand('~/R/cf_model_comparison_20261013_R46'),path.expand('~/R/x86_64-pc-linux-gnu-library/4.5'),.libPaths()))
library(glmmTMB)
out<-'age_sensitivity_20261119'
comparisons<-fixed<-list()
for(trait in c('clutch','fledging')) for(kind in c('age_class','exact_age_linear','exact_age_quadratic')) {
 path<-file.path(out,paste0(trait,'_',kind,'.rds'))
 while(!file.exists(path))Sys.sleep(15)
 a<-readRDS(path)
 if(!isTRUE(a$sdr$pdHess))next
 d<-model.frame(a)
 if(trait=='fledging') {
   y<-model.response(d);d$fledged<-y[,1];d$clutch<-rowSums(y)
 }
 cf<-fixef(a)
 start<-a$obj$env$parList(par=a$fit$parfull)
 start<-start[intersect(c('beta','betazi','betadisp','theta','thetazi','thetadisp'),names(start))]
 start<-start[lengths(start)>0]
 V<-a$sdr$cov.fixed
 optimizer<-function(par,fn,gr,control=list(),...) {
   eig<-eigen((V+t(V))/2,symmetric=TRUE)
   stopifnot(min(eig$values)>0,length(par)==nrow(V))
   S<-eig$vectors%*%diag(sqrt(eig$values))
   ff<-function(q)fn(drop(par+S%*%q))
   gg<-function(q)drop(crossprod(S,gr(drop(par+S%*%q))))
   if(kind=='age_class')z<-optim(rep(0,length(par)),ff,gg,method='BFGS',control=list(maxit=1500,reltol=1e-12)) else {
     z<-nlminb(rep(0,length(par)),ff,gg,control=list(iter.max=1500,eval.max=3000,rel.tol=1e-11))
     z$value<-z$objective
   }
   list(par=drop(par+S%*%z$par),objective=z$value,convergence=z$convergence,message=z$message)
 }
 refined<-file.path(out,paste0(trait,'_',kind,'_refined.rds'))
 if(file.exists(refined))b<-readRDS(refined) else {
   b<-glmmTMB(formula(a),dispformula=formula(a,component='disp'),ziformula=formula(a,component='zi'),
     family=family(a),data=d,start=start,control=glmmTMBControl(optimizer=optimizer,parallel=4))
   saveRDS(b,refined)
 }
 gradient<-max(abs(b$sdr$gradient.fixed*sqrt(diag(b$sdr$cov.fixed))))
 dg<-data.frame(trait,model=kind,n=nobs(b),convergence=b$fit$convergence,pdHess=isTRUE(b$sdr$pdHess),
    logLik=as.numeric(logLik(b)),AIC=AIC(b),scaled_gradient=gradient,logLik_difference=as.numeric(logLik(b)-logLik(a)))
 comparisons[[paste(trait,kind)]]<-dg
 write.csv(do.call(rbind,comparisons),file.path(out,'validated_diagnostics.csv'),row.names=FALSE)
 print(dg);flush.console()
 if(b$fit$convergence==0 && isTRUE(b$sdr$pdHess) && gradient<.005 && abs(dg$logLik_difference)<.01) {
  co<-summary(b)$coefficients
  fixed[[paste(trait,kind)]]<-do.call(rbind,lapply(names(co),function(component){x<-co[[component]];if(!length(x))return(NULL)
    data.frame(trait,model=kind,component,term=rownames(x),estimate=x[,1],SE=x[,2],lower=x[,1]-qnorm(.975)*x[,2],upper=x[,1]+qnorm(.975)*x[,2])}))
  write.csv(do.call(rbind,fixed),file.path(out,'validated_fixed_effects.csv'),row.names=FALSE)
 }
}
cat('AGE VALIDATION COMPLETE\n')
