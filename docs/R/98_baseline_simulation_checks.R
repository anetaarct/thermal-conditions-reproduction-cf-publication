# Conditional plug-in checks for regenerated baseline models; no held-out validation.
library(glmmTMB)
base<-Sys.getenv('WB_MODEL_DIR','models/within_between_20261026')
checks<-list();B<-300L
pmf_compois<-function(mu,phi) {
 nu<-1/phi;k<-0:max(40,ceiling(max(mu)+20));lf<-lgamma(k+1)
 probabilities<-function(a){z<-outer(a,k)-outer(nu,lf);z<-exp(z-do.call(pmax,as.data.frame(z)));z/rowSums(z)}
 lo<-nu*log(pmax(mu-3,.001))-20;hi<-nu*log(mu+3)+20
 for(i in 1:65){mid<-(lo+hi)/2;p<-probabilities(mid);low<-drop(p%*%k)<mu;lo[low]<-mid[low];hi[!low]<-mid[!low]}
 p<-probabilities((lo+hi)/2)
 stopifnot(max(abs(drop(p%*%k)-mu))<1e-7,max(p[,ncol(p)])<1e-10)
 p
}
for(trait in c('clutch','fledging')) {
 set.seed(20261104+match(trait,c('clutch','fledging')))
 m<-readRDS(file.path(base,paste0(trait,'_baseline_selected.rds')))
 d<-readRDS(file.path(base,paste0(trait,'_analysis_data.rds')))
 mu<-predict(m,type=if(trait=='clutch')'response' else 'conditional');phi<-predict(m,type='disp')
 if(trait=='clutch') {
  p<-pmf_compois(mu,phi);cdf<-t(apply(p,1,cumsum))
  y<-replicate(B,rowSums(cdf<runif(nrow(d))))
  observed<-d$clutch
 } else {
  zi<-predict(m,type='zprob')
  y<-replicate(B,{z<-rbinom(nrow(d),d$clutch,rbeta(nrow(d),mu*phi,(1-mu)*phi));z[runif(nrow(d))<zi]<-0;z})
  y<-y/d$clutch;observed<-d$fledged/d$clutch
 }
 groups<-c(list(overall=seq_len(nrow(d))),split(seq_len(nrow(d)),d$year_f),split(seq_len(nrow(d)),d$age_class))
 for(g in names(groups)) {
  ix<-groups[[g]];z<-y[ix,,drop=FALSE];obs<-observed[ix]
  stats<-list(mean=colMeans(z),variance=apply(z,2,var));v<-c(mean=mean(obs),variance=var(obs))
  if(trait=='fledging') {stats$complete_failure<-colMeans(z==0);stats$complete_success<-colMeans(z==1);v<-c(v,complete_failure=mean(obs==0),complete_success=mean(obs==1))}
  for(name in names(stats)) checks[[length(checks)+1L]]<-data.frame(trait,group=g,statistic=name,
   observed=unname(v[name]),predicted=mean(stats[[name]]),lower=quantile(stats[[name]],.025),upper=quantile(stats[[name]],.975),simulations=B)
 }
}
dir.create('tables/reproduced',recursive=TRUE,showWarnings=FALSE)
write.csv(do.call(rbind,checks),'tables/reproduced/baseline_simulation_checks.csv',row.names=FALSE)
