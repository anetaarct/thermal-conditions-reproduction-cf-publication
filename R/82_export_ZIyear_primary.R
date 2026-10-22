.libPaths(c(path.expand('~/R/cf_model_comparison_20261013_R46'),path.expand('~/R/x86_64-pc-linux-gnu-library/4.5'),.libPaths()))
suppressPackageStartupMessages({library(glmmTMB);library(MASS)})
setwd('/home/arct/CF_model_comparison_20261013');out<-'zi_year_primary_20261021';dir.create(out,showWarnings=FALSE)
p<-'zi_year_sensitivity_20261020/zi_year_selected.rds';h<-tools::md5sum(p);m<-readRDS(p);d<-readRDS('glmmTMB_comparison_20261016/analysis_data.rds');stopifnot(nobs(m)==9790,m$fit$convergence==0,isTRUE(m$sdr$pdHess))
put<-function(x,n)write.csv(x,file.path(out,n),row.names=FALSE)
b<-fixef(m);bn<-c(names(b$cond),paste0('zi~',names(b$zi)),paste0('disp~',names(b$disp)));bb<-c(unname(b$cond),unname(b$zi),unname(b$disp));V<-vcov(m,full=TRUE)[bn,bn,drop=FALSE];stopifnot(min(eigen(V,symmetric=TRUE,only.values=TRUE)$values)>0);set.seed(20261021);dr<-MASS::mvrnorm(8000,bb,V);nc<-length(b$cond);nz<-length(b$zi)
ref<-function(n){nd<-d[rep(1,n),,drop=FALSE];nd$rel_LD<-0;nd$temp_mean_c<-0;nd$year_c<-0;nd$age_class<-factor('YOUNG',levels=levels(d$age_class));nd}
grids<-list();a<-ref(100);a$x<-seq(quantile(d$temp_mean_c,.05),quantile(d$temp_mean_c,.95),length.out=100);a$temp_mean_c<-a$x;grids$temperature<-a
a<-ref(200);a$x<-rep(seq(quantile(d$rel_LD,.05),quantile(d$rel_LD,.95),length.out=100),2);a$rel_LD<-a$x;a$age_class<-factor(rep(c('YOUNG','OLD'),each=100),levels=levels(d$age_class));grids$date<-a
a<-ref(39);a$x<-1981:2019;a$year_c<-a$x-unique(as.numeric(as.character(d$year_f))-d$year_c);grids$year<-a
a<-ref(2);a$x<-1:2;a$age_class<-factor(c('YOUNG','OLD'),levels=levels(d$age_class));grids$age<-a
for(kind in names(grids)){
 nd<-grids[[kind]];X<-model.matrix(~rel_LD+age_class+temp_mean_c+year_c,nd);mu<-plogis(drop(X%*%b$cond));phi<-exp(drop(X%*%b$disp));zi<-plogis(b$zi[1]);success<-(1-zi)*mu
 md<-plogis(dr[,1:nc,drop=FALSE]%*%t(X));pd<-exp(dr[,nc+nz+seq_len(length(b$disp)),drop=FALSE]%*%t(X));zd<-plogis(dr[,nc+1]);sd<-(1-zd)*md
 # Analytic fixed designs set all random effects to zero, including dispersion and zero components.
 zero_nd<-nd;zero_nd$year_f<-NA;zero_nd$female_id<-NA;stopifnot(max(abs(predict(m,newdata=zero_nd,type='response',re.form=NA)-success))<1e-7,max(abs(predict(m,newdata=zero_nd,type='zprob',re.form=NA)-zi))<1e-7,max(abs(predict(m,newdata=zero_nd,type='disp',re.form=NA)-phi))<1e-7)
 for(component in c('success','phi')){val<-if(component=='success')sd else pd;point<-if(component=='success')success else phi;qq<-t(apply(val,2,quantile,probs=c(.025,.975)));put(data.frame(x=nd$x,age_class=nd$age_class,estimate=point,CI_low=qq[,1],CI_high=qq[,2]),paste0('ZIyear_',kind,'_',component,'.csv'))}
}
co<-summary(m)$coefficients;fx<-do.call(rbind,lapply(names(co),function(cmp){a<-co[[cmp]];data.frame(component=cmp,term=rownames(a),estimate=a[,1],SE=a[,2],lower=a[,1]-qnorm(.975)*a[,2],upper=a[,1]+qnorm(.975)*a[,2])}));put(fx,'ZIyear_fixed.csv')
# Log-SD Wald intervals, including fixed/random-SD covariance from the fitted model Hessian.
par<-m$fit$par;se<-sqrt(diag(m$sdr$cov.fixed));ix<-which(names(par)%in%c('theta','thetazi','thetadisp'));print(data.frame(name=names(par)[ix],value=par[ix],SE=se[ix]));re<-VarCorr(m);rr<-do.call(rbind,lapply(names(re),function(cmp)do.call(rbind,lapply(names(re[[cmp]]),function(gr)data.frame(component=cmp,group=gr,SD=attr(re[[cmp]][[gr]],'stddev'))))))
# Match components explicitly to internal parameter vector order.
parts<-list(cond=which(names(par)=='theta'),zi=which(names(par)=='thetazi'),disp=which(names(par)=='thetadisp'))
rs<-do.call(rbind,lapply(names(parts),function(cmp){a<-rr[rr$component==cmp,,drop=FALSE];ii<-parts[[cmp]];stopifnot(length(ii)==nrow(a),max(abs(exp(par[ii])-a$SD))<1e-7);a$lower<-exp(par[ii]-qnorm(.975)*se[ii]);a$upper<-exp(par[ii]+qnorm(.975)*se[ii]);a}));put(rs,'ZIyear_random_SD_intervals.csv')
z<-b$zi[1];zs<-sqrt(vcov(m)$zi[1,1]);put(data.frame(target='extra-zero probability at zero year effect; not annual average or total failure probability',estimate=plogis(z),lower=plogis(z-qnorm(.975)*zs),upper=plogis(z+qnorm(.975)*zs)),'ZIyear_reference_zi.csv')
put(data.frame(model='ZI_year_primary_ML',n=nobs(m),logLik=as.numeric(logLik(m)),AIC=AIC(m),convergence=m$fit$convergence,pdHess=m$sdr$pdHess,scaled_gradient=max(abs(m$sdr$gradient.fixed*sqrt(diag(m$sdr$cov.fixed))))),'ZIyear_diagnostics.csv')
stopifnot(identical(h,tools::md5sum(p)));capture.output(sessionInfo(),file=file.path(out,'session_info.txt'));writeLines(c('8000 joint asymptotic coefficient draws; point estimates are ML, not draw medians.','All random effects zero for curves, including zero-year and precision-year. Temperature x centred at 5.098595980711215 C.','Original selected fit unchanged.'),file.path(out,'COMPLETE.txt'));print(rs);cat('COMPLETE\n')
