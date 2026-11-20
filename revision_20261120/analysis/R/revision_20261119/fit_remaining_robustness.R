.libPaths(c(path.expand('~/R/cf_model_comparison_20261013_R46'),path.expand('~/R/x86_64-pc-linux-gnu-library/4.5'),.libPaths()))
Sys.setenv(OMP_NUM_THREADS=1,OPENBLAS_NUM_THREADS=1)
suppressPackageStartupMessages({library(glmmTMB);library(MASS)})
task<-Sys.getenv('ROBUSTNESS_TASK','all')
out<-file.path('remaining_robustness_20261120',task);dir.create(out,recursive=TRUE,showWarnings=FALSE)
contrast_script<-if(file.exists('export_sensitivity_contrasts.R'))'export_sensitivity_contrasts.R' else 'R/revision_20261119/export_sensitivity_contrasts.R'
code<-readLines(contrast_script);eval(parse(text=code[seq_len(grep('^age<-list',code)[1]-1)]))
diagrows<-effects<-samples<-fixed<-random<-comparisons<-list()
flush<-function(){for(n in c('diagrows','effects','samples','fixed','random','comparisons')){x<-get(n);if(length(x))write.csv(do.call(rbind,x),file.path(out,paste0(n,'.csv')),row.names=FALSE)}}
key<-function(x)paste(x,ave(seq_along(x),x,FUN=seq_along),sep='_')
for(check in c('repeated_females','scale_year'))for(trait in c('clutch','fledging'))for(kind in c('baseline','within_between')){
 if(check=='scale_year'&&trait=='fledging')next # Already present in the primary fledging model.
 tag<-paste(check,trait,kind,sep='_');if(task!='all'&&task!=tag)next
 if(tag=='repeated_females_clutch_within_between'&&Sys.getenv('RUN_UNVALIDATED_CLUTCH_WB','false')!='true'){
  writeLines('This attempted specification did not yield a numerically validated fit in the revision. No inference is reported. See the archived objective and optimiser audits. Set RUN_UNVALIDATED_CLUTCH_WB=true only to rerun the unvalidated attempt.',file.path(out,'clutch_within_between_status.txt'))
  next
 }
 cat(format(Sys.time()),tag,'START\n');flush.console()
 reference<-readRDS(file.path('models/clean_20261119',paste0(trait,'_',kind,'_selected.rds')))
 d<-readRDS(file.path('models/clean_20261119',paste0(trait,'_analysis_data.rds')))
 if(check=='repeated_females'){
  counts<-table(d$female_id);d<-droplevels(d[as.character(d$female_id)%in%names(counts[counts>1]),])
  # Preserve full-sample centring/decomposition so the sensitivity changes only inclusion.
  stopifnot(all(table(d$female_id)>=2))
 }
 samples[[tag]]<-data.frame(check,trait,model=kind,n=nrow(d),females=length(unique(d$female_id)),seasons=length(unique(d$year)))
 fits<-list()
 methods<-c('NLMINB','BFGS')
 for(method in methods){
  seed<-if(length(fits))fits[[length(fits)]] else reference
  lifted<-check=='repeated_females'&&trait=='clutch'&&kind=='within_between'&&!length(fits)
  if(lifted)seed<-readRDS(file.path(if(task=='all')out else 'remaining_robustness_20261120/repeated_females_clutch_baseline','repeated_females_clutch_baseline_selected.rds'))
  p<-seed$obj$env$parList(par=seed$fit$parfull)
  start<-p[intersect(c('beta','betazi','betadisp','theta','thetazi','thetadisp'),names(p))];start<-start[lengths(start)>0]
  modes<-ranef(seed)
  start$b<-c(modes$cond$female_id[as.character(levels(d$female_id)),1],modes$cond$year_f[as.character(levels(d$year_f)),1])
  if(length(p$bdisp)&&!is.null(modes$disp$year_f))start$bdisp<-modes$disp$year_f[as.character(levels(d$year_f)),1]
  if(length(p$bzi)&&!is.null(modes$zi$year_f))start$bzi<-modes$zi$year_f[as.character(levels(d$year_f)),1]
  stopifnot(all(is.finite(start$b)))
  if(check=='scale_year'&&!length(start$thetadisp))start$thetadisp<-log(.05)
  V<-seed$sdr$cov.fixed;oldkeys<-key(rownames(V))
  if(lifted){
   columns<-colnames(model.matrix(~rel_LD+temp_mean_c+year_c+age_class,d))
   active<-colnames(model.matrix(~rel_LD+temp_within+temp_between+year_c+age_class,d))
   cf<-fixef(seed)
   for(component in c('cond','disp')){
    b<-cf[[component]];b['temp_within']<-b['temp_mean_c'];b['temp_between']<-b['temp_mean_c']
    start[[if(component=='cond')'beta' else 'betadisp']]<-unname(b[active])
   }
   oldnames<-rownames(V);newnames<-c(rep('beta',length(active)),rep('betadisp',length(active)),rep('theta',2))
   A<-matrix(0,length(newnames),length(oldnames));D<-matrix(0,length(newnames),2)
   for(g in c('beta','betadisp','theta')){
    ni<-which(newnames==g);oi<-which(oldnames==g)
    if(g=='theta')A[cbind(ni,oi)]<-1 else{
     mapped<-match(sub('temp_(within|between)','temp_mean_c',active),columns)
     A[cbind(ni,oi[mapped])]<-1
     D[ni[active=='temp_within'],if(g=='beta')1 else 2]<-.5
     D[ni[active=='temp_between'],if(g=='beta')1 else 2]<--.5
    }
   }
   V<-A%*%V%*%t(A)+D%*%diag(c(.01^2,.05^2))%*%t(D)
   rownames(V)<-newnames;oldkeys<-key(newnames)
  }
  stopifnot(length(oldkeys)==nrow(V))
  optimizer<-function(par,fn,gr,control=list(),...){
   if(trait=='clutch'){environment(fn)$inner.control$maxit<-100L;environment(fn)$inner.control$grad.tol<-1e-7}
   newkeys<-key(names(par));idx<-match(newkeys,oldkeys)
   stopifnot(sum(!is.na(idx))>=length(par)-1)
   C<-diag(.25,length(par));ok<-which(!is.na(idx));C[ok,ok]<-V[idx[ok],idx[ok],drop=FALSE]
   C[!is.finite(C)]<-0
   e<-eigen((C+t(C))/2,symmetric=TRUE);S<-e$vectors%*%diag(sqrt(pmax(e$values,1e-8)))
   ff<-function(q)fn(drop(par+S%*%q));gg<-function(q)drop(crossprod(S,gr(drop(par+S%*%q))))
   if(method=='BFGS')z<-optim(rep(0,length(par)),ff,gg,method='BFGS',control=list(maxit=2500,reltol=1e-12,trace=1,REPORT=20)) else{
    z<-nlminb(rep(0,length(par)),ff,gg,control=list(iter.max=2500,eval.max=5000,rel.tol=1e-11,step.max=.25,trace=1));z$value<-z$objective
   }
   list(par=drop(par+S%*%z$par),objective=z$value,convergence=z$convergence,message=z$message)
  }
  file<-file.path(out,paste0(tag,'_',method,'.rds'))
  m<-tryCatch({if(file.exists(file))readRDS(file) else{
   if(check=='scale_year')ans<-update(seed,data=d,dispformula=update(formula(reference,component='disp'),~.+(1|year_f)),start=start,control=glmmTMBControl(optimizer=optimizer,parallel=1)) else ans<-update(seed,formula=formula(reference),dispformula=formula(reference,component='disp'),data=d,start=start,control=glmmTMBControl(optimizer=optimizer,parallel=1))
   saveRDS(ans,file);ans}},error=function(e){cat('ERROR',tag,method,conditionMessage(e),'\n');NULL})
  if(is.null(m))next
  dg<-data.frame(check,trait,model=kind,optimizer=method,n=nobs(m),convergence=m$fit$convergence,pdHess=isTRUE(m$sdr$pdHess),logLik=as.numeric(logLik(m)),objective=m$fit$objective,AIC=AIC(m),scaled_gradient=tryCatch(max(abs(m$sdr$gradient.fixed*sqrt(diag(m$sdr$cov.fixed)))),error=function(e)NA_real_))
  diagrows[[paste(tag,method)]]<-dg;fits[[method]]<-m;flush();print(dg);flush.console()
 }
 valid<-vapply(fits,function(m)isTRUE(m$sdr$pdHess)&&m$fit$convergence==0&&is.finite(as.numeric(logLik(m))),logical(1))
 if(!any(valid)){cat(tag,'NO VALID FIT\n');next}
 candidates<-fits[valid];m<-candidates[[which.max(vapply(candidates,function(x)as.numeric(logLik(x)),numeric(1)))]]
 saveRDS(m,file.path(out,paste0(tag,'_selected.rds')))
 predictors<-if(kind=='baseline')'temp_mean_c' else c('temp_within','temp_between')
 effects[[tag]]<-do.call(rbind,lapply(predictors,function(p)cbind(data.frame(check,trait,model=kind,predictor=p),evaluate_contrast(m,trait,p,20261120)[,c('mean_change','lower','upper','consistency_ratio','ratio_lower','ratio_upper')])))
 co<-summary(m)$coefficients
 fixed[[tag]]<-do.call(rbind,lapply(names(co),function(component){x<-co[[component]];if(!length(x))return(NULL);data.frame(check,trait,model=kind,component,term=rownames(x),estimate=x[,1],SE=x[,2],lower=x[,1]-1.96*x[,2],upper=x[,1]+1.96*x[,2])}))
 re<-VarCorr(m);random[[tag]]<-do.call(rbind,lapply(names(re),function(component)do.call(rbind,lapply(names(re[[component]]),function(group)data.frame(check,trait,model=kind,component,group,SD=unname(attr(re[[component]][[group]],'stddev')))))))
 comparisons[[tag]]<-data.frame(check,trait,model=kind,reference_n=nobs(reference),sensitivity_n=nobs(m),reference_AIC=AIC(reference),sensitivity_AIC=AIC(m),AIC_comparable=check=='scale_year',delta_AIC=if(check=='scale_year')AIC(m)-AIC(reference) else NA_real_,optimizer_logLik_range=diff(range(vapply(candidates,function(x)as.numeric(logLik(x)),numeric(1)))))
 capture.output(summary(m),file=file.path(out,paste0(tag,'_summary.txt')))
 flush();cat(tag,'COMPLETE\n');flush.console()
}
capture.output(sessionInfo(),file=file.path(out,'session_info.txt'))
cat('REMAINING ROBUSTNESS COMPLETE\n')
