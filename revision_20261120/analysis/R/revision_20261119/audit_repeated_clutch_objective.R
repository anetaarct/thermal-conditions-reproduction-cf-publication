.libPaths(c(path.expand('~/R/cf_model_comparison_20261013_R46'),path.expand('~/R/x86_64-pc-linux-gnu-library/4.5'),.libPaths()))
library(glmmTMB)
m<-readRDS('remaining_robustness_20261120/repeated_females_clutch_baseline/repeated_females_clutch_baseline_selected.rds')
d<-model.frame(m)
raw<-readRDS('models/clean_20261119/clutch_analysis_data.rds')
key<-paste(raw$female_id,raw$year);idx<-match(paste(d$female_id,d$year_f),key)
stopifnot(!anyNA(idx));d$temp_within<-raw$temp_within[idx];d$temp_between<-raw$temp_between[idx]
reference<-readRDS('models/clean_20261119/clutch_within_between_selected.rds')
cols<-colnames(model.matrix(~rel_LD+temp_within+temp_between+year_c+age_class,d))
p<-m$obj$env$parList(par=m$fit$parfull);st<-p[intersect(c('beta','betadisp','theta','b'),names(p))]
for(component in c('cond','disp')){b<-fixef(m)[[component]];b['temp_within']<-b['temp_mean_c'];b['temp_between']<-b['temp_mean_c'];st[[if(component=='cond')'beta' else 'betadisp']]<-unname(b[cols])}
optimizer<-function(par,fn,gr,control=list(),...){
 inner<-environment(fn)$inner.control;inner$maxit<-1000L;inner$grad.tol<-1e-9;environment(fn)$inner.control<-inner
 V<-m$sdr$cov.fixed;oldnames<-rownames(V);newnames<-names(par)
 basecols<-colnames(model.matrix(~rel_LD+temp_mean_c+year_c+age_class,d))
 A<-matrix(0,length(par),length(oldnames));D<-matrix(0,length(par),2)
 for(g in c('beta','betadisp','theta')){ni<-which(newnames==g);oi<-which(oldnames==g);if(g=='theta')A[cbind(ni,oi)]<-1 else{mapped<-match(sub('temp_(within|between)','temp_mean_c',cols),basecols);A[cbind(ni,oi[mapped])]<-1;D[ni[cols=='temp_within'],if(g=='beta')1 else 2]<-.5;D[ni[cols=='temp_between'],if(g=='beta')1 else 2]<--.5}}
 V<-A%*%V%*%t(A)+D%*%diag(c(.01^2,.05^2))%*%t(D)
 e<-eigen(V,symmetric=TRUE);S<-e$vectors%*%diag(sqrt(pmax(e$values,1e-8)))
 q<-c(.0792275,.464835,-.425393,-1.48152,.0786817,.394192,-.146614,.120892,.0251549,-.286053,-.00889812,.328108,.172635,-.00158995)
 par<-drop(par+S%*%q)
 cat('Target parameters',par,'\n');flush.console()
 for(i in 1:3){t<-system.time(v<-fn(par));g<-gr(par);cat(i,'objective',v,'maxgradient',max(abs(g)),'elapsed',t[3],'\n');flush.console()}
 stop('Objective audit complete; no fit requested')
}
try(update(m,formula=formula(reference),dispformula=formula(reference,component='disp'),family=compois(),data=d,start=st,control=glmmTMBControl(optimizer=optimizer,parallel=1)))
