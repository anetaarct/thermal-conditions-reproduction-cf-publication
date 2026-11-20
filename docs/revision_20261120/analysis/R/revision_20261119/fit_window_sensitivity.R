.libPaths(c(path.expand('~/R/cf_model_comparison_20261013_R46'),path.expand('~/R/x86_64-pc-linux-gnu-library/4.5'), .libPaths()))
Sys.setenv(OMP_NUM_THREADS=1,OPENBLAS_NUM_THREADS=1)
library(glmmTMB)
library(MASS)
out <- 'window_sensitivity_20261119'
dir.create(out,showWarnings=FALSE)
weather <- read.csv('data_derived/smhi_hoburg_daily_temperature.csv')
weather$date <- as.Date(weather$date)
windows <- data.frame(window=c('original','apr01_may07','apr01_apr30','apr15_may07'),
  start=c('03-26','04-01','04-01','04-15'),end=c('05-07','05-07','04-30','05-07'))
annual <- do.call(rbind,lapply(seq_len(nrow(windows)),function(i)do.call(rbind,lapply(1981:2025,function(y){
  a<-as.Date(paste(y,windows$start[i],sep='-'));b<-as.Date(paste(y,windows$end[i],sep='-'))
  dates<-seq(a,b,by='day');x<-weather$temperature[match(dates,weather$date)]
  stopifnot(length(x)==as.integer(b-a)+1L,!anyNA(x))
  data.frame(window=windows$window[i],year=y,start=a,end=b,days=length(x),temperature=mean(x))
}))))
write.csv(annual,file.path(out,'annual_exposures.csv'),row.names=FALSE)
# Reuse only the contrast function; do not execute unrelated export loops.
contrast_script<-if(file.exists('export_sensitivity_contrasts.R'))'export_sensitivity_contrasts.R' else 'R/revision_20261119/export_sensitivity_contrasts.R'
code<-readLines(contrast_script)
cut<-grep('^age<-list',code)[1]
eval(parse(text=code[seq_len(cut-1)]))
all_diag<-all_fixed<-all_contrasts<-list()
flush_tables<-function(){
  if(length(all_diag))write.csv(do.call(rbind,all_diag),file.path(out,'diagnostics.csv'),row.names=FALSE)
  if(length(all_fixed))write.csv(do.call(rbind,all_fixed),file.path(out,'fixed_effects.csv'),row.names=FALSE)
  if(length(all_contrasts))write.csv(do.call(rbind,all_contrasts),file.path(out,'reference_contrasts.csv'),row.names=FALSE)
}
for(trait in c('clutch','fledging'))for(kind in c('baseline','within_between')){
  reference<-readRDS(file.path('models/clean_20261119',paste0(trait,'_',kind,'_selected.rds')))
  V<-reference$sdr$cov.fixed
  E<-eigen((V+t(V))/2,symmetric=TRUE);stopifnot(min(E$values)>0)
  S<-E$vectors%*%diag(sqrt(E$values))
  make_optimizer<-function(method)function(par,fn,gr,control=list(),...){
    stopifnot(length(par)==nrow(S))
    if(trait=='clutch'){
      environment(fn)$inner.control$maxit<-100L
      environment(fn)$inner.control$grad.tol<-1e-7
    }
    ff<-function(q)fn(drop(par+S%*%q));gg<-function(q)drop(crossprod(S,gr(drop(par+S%*%q))))
    if(method!='BFGS')stopifnot(abs(ff(rep(0,length(par)))+as.numeric(logLik(seed)))<.01)
    if(method=='BFGS'){
      z<-optim(rep(0,length(par)),ff,gg,method='BFGS',control=list(maxit=1500,reltol=1e-12))
      list(par=drop(par+S%*%z$par),objective=z$value,convergence=z$convergence,message=z$message)
    }else{
      z<-nlminb(rep(0,length(par)),ff,gg,control=list(iter.max=1500,eval.max=3000,rel.tol=1e-10))
      list(par=drop(par+S%*%z$par),objective=z$objective,convergence=z$convergence,message=z$message)
    }
  }
  for(w in windows$window){
    tag<-paste(trait,kind,w,sep='_');cat(format(Sys.time()),tag,'START\n');flush.console()
    d<-readRDS(file.path('models/clean_20261119',paste0(trait,'_analysis_data.rds')))
    a<-annual[annual$window==w,];exposure<-a$temperature[match(d$year,a$year)]
    stopifnot(!anyNA(exposure))
    if(w=='original')stopifnot(max(abs(exposure-d$temp_mean))<1e-8)
    d$temp_mean<-exposure;d$temp_mean_c<-exposure-mean(exposure)
    d$temp_between<-ave(d$temp_mean_c,d$female_id,FUN=mean)
    d$temp_within<-d$temp_mean_c-d$temp_between
    saveRDS(d,file.path(out,paste0(tag,'_data.rds')))
    write.csv(d,file.path(out,paste0(tag,'_data.csv')),row.names=FALSE)
    fits<-list()
    methods<-c('BFGS','NLMINB')
    for(method in methods){
      file<-file.path(out,paste0(tag,'_',method,'.rds'))
      if(w!='original'&&file.exists(file)){
        candidate<-readRDS(file)
        if(!isTRUE(candidate$sdr$pdHess)||candidate$fit$convergence!=0){
          archived<-sub('[.]rds$',paste0('_inadmissible_',format(Sys.time(),'%Y%m%d_%H%M%S'),'.rds'),file)
          stopifnot(file.rename(file,archived))
        }
      }
      if(w=='original')m<-reference else if(file.exists(file))m<-readRDS(file) else{
        seed<-if(length(fits))fits[[1]] else reference
        if(method!='BFGS'){
          seed_V<-seed$sdr$cov.fixed
          seed_E<-eigen((seed_V+t(seed_V))/2,symmetric=TRUE)
          stopifnot(min(seed_E$values)>0)
          S<-seed_E$vectors%*%diag(sqrt(seed_E$values))
        }
        p<-seed$obj$env$parList(par=seed$fit$parfull)
        start_names<-c('beta','betazi','betadisp','theta','thetazi','thetadisp')
        if(length(fits))start_names<-c(start_names,'b','bzi','bdisp')
        start<-p[intersect(names(p),start_names)]
        start<-start[lengths(start)>0]
        rhs<-paste('rel_LD +',if(kind=='baseline')'temp_mean_c' else 'temp_within + temp_between','+ year_c + age_class')
        m<-glmmTMB(as.formula(paste(if(trait=='clutch')'clutch' else 'cbind(fledged,clutch-fledged)','~',rhs,'+ (1|female_id)+(1|year_f)')),
          dispformula=as.formula(paste('~',rhs,if(trait=='fledging')'+(1|year_f)' else '')),
          ziformula=if(trait=='fledging')~1+(1|year_f) else ~0,
          family=if(trait=='clutch')compois() else betabinomial(),data=d,start=start,
          control=glmmTMBControl(optimizer=make_optimizer(method),parallel=4))
        saveRDS(m,file)
      }
      dg<-data.frame(trait,model=kind,window=w,optimizer=if(w=='original')'archived_validated' else method,
        n=nobs(m),convergence=m$fit$convergence,pdHess=isTRUE(m$sdr$pdHess),logLik=as.numeric(logLik(m)),AIC=AIC(m),
        scaled_gradient=max(abs(m$sdr$gradient.fixed*sqrt(diag(m$sdr$cov.fixed)))))
      all_diag[[paste(tag,method)]]<-dg;flush_tables();print(dg);flush.console()
      if(!dg$pdHess||dg$convergence!=0||!is.finite(dg$scaled_gradient)||dg$scaled_gradient>.005)stop('Invalid fit: ',tag,method)
      fits[[method]]<-m
      if(w=='original')break
    }
    if(length(fits)==2)stopifnot(abs(as.numeric(logLik(fits[[1]]))-as.numeric(logLik(fits[[2]])))<.001)
    m<-fits[[which.max(vapply(fits,function(x)as.numeric(logLik(x)),numeric(1)))]]
    saveRDS(m,file.path(out,paste0(tag,'_selected.rds')))
    co<-summary(m)$coefficients
    all_fixed[[tag]]<-do.call(rbind,lapply(names(co),function(component){x<-co[[component]];if(!length(x))return(NULL)
      data.frame(trait,model=kind,window=w,component,term=rownames(x),estimate=x[,1],SE=x[,2],lower=x[,1]-1.96*x[,2],upper=x[,1]+1.96*x[,2])}))
    predictors<-if(kind=='baseline')'temp_mean_c' else c('temp_within','temp_between')
    all_contrasts[[tag]]<-do.call(rbind,lapply(predictors,function(p)cbind(data.frame(trait,model=kind,window=w,predictor=p),
      evaluate_contrast(m,trait,p,20261122)[,c('mean_change','lower','upper','consistency_ratio','ratio_lower','ratio_upper')])) )
    flush_tables();cat(tag,'COMPLETE\n');flush.console()
  }
}
capture.output(sessionInfo(),file=file.path(out,'session_info.txt'))
cat('WINDOW SENSITIVITY COMPLETE\n')
