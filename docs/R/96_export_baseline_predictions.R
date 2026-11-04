# Portable baseline temperature/date/age exports from newly fitted public-data models.
suppressPackageStartupMessages({library(glmmTMB);library(MASS)})
base <- Sys.getenv('WB_MODEL_DIR','models/within_between_20261026')
out <- 'tables/reproduced/baseline_predictions'
dir.create(out,recursive=TRUE,showWarnings=FALSE)
curves <- contrasts <- list()
for(trait in c('clutch','fledging')) {
  m <- readRDS(file.path(base,paste0(trait,'_baseline_selected.rds')))
  d <- readRDS(file.path(base,paste0(trait,'_analysis_data.rds')))
  stopifnot(m$fit$convergence==0,isTRUE(m$sdr$pdHess))
  cf <- fixef(m); nc<-length(cf$cond);nd<-length(cf$disp)
  names_full<-c(names(cf$cond),paste0('disp~',names(cf$disp)),if(length(cf$zi))paste0('zi~',names(cf$zi)))
  V<-vcov(m,full=TRUE)[names_full,names_full,drop=FALSE]
  set.seed(20261104+match(trait,c('clutch','fledging')))
  draws<-MASS::mvrnorm(4000,c(unname(cf$cond),unname(cf$disp),unname(cf$zi)),V)
  evaluate <- function(g) {
    formula <- if(trait=='clutch') ~rel_LD+temp_mean_c+year_c+age_class else ~rel_LD+age_class+temp_mean_c+year_c
    X<-model.matrix(formula,g)
    eta<-drop(X%*%cf$cond); etad<-X%*%t(draws[,seq_len(nc),drop=FALSE])
    phi<-exp(drop(X%*%cf$disp)); phid<-exp(X%*%t(draws[,nc+seq_len(nd),drop=FALSE]))
    if(trait=='clutch') {
      mu<-exp(eta); mud<-exp(etad); consistency<-sqrt(compois()$variance(mu,phi))
      cd<-matrix(sqrt(compois()$variance(as.vector(mud),as.vector(phid))),nrow(g),nrow(draws))
    } else {
      mu<-(1-plogis(cf$zi[1]))*plogis(eta)
      mud<-sweep(plogis(etad),2,1-plogis(draws[,nc+nd+1]),'*')
      consistency<-phi;cd<-phid
    }
    list(mean=mu,means=mud,consistency=consistency,consistencies=cd)
  }
  grid <- function(n) data.frame(rel_LD=rep(0,n),temp_mean_c=0,year_c=0,
    age_class=factor(rep('YOUNG',n),levels=c('YOUNG','OLD')))
  for(predictor in c('temp_mean_c','rel_LD','age_class')) {
    g<-grid(if(predictor=='age_class')2 else 101)
    g[[predictor]]<-if(predictor=='age_class')factor(c('YOUNG','OLD'),levels=c('YOUNG','OLD')) else
      seq(quantile(d[[predictor]],.05),quantile(d[[predictor]],.95),length.out=101)
    values<-evaluate(g)
    for(component in c('mean','consistency')) {
      v<-values[[if(component=='mean')'means' else 'consistencies']]
      ci<-t(apply(v,1,quantile,c(.025,.975)))
      curves[[paste(trait,predictor,component)]]<-data.frame(trait,predictor,component,
        x=as.character(g[[predictor]]),estimate=values[[component]],lower=ci[,1],upper=ci[,2])
    }
    g<-grid(2)
    g[[predictor]]<-if(predictor=='age_class')factor(c('YOUNG','OLD'),levels=c('YOUNG','OLD')) else c(0,1)
    v<-evaluate(g);delta<-v$means[2,]-v$means[1,];ratio<-v$consistencies[2,]/v$consistencies[1,]
    contrasts[[paste(trait,predictor)]]<-data.frame(trait,predictor,
      mean_change=v$mean[2]-v$mean[1],mean_lower=quantile(delta,.025),mean_upper=quantile(delta,.975),
      consistency_ratio=v$consistency[2]/v$consistency[1],ratio_lower=quantile(ratio,.025),ratio_upper=quantile(ratio,.975))
  }
}
write.csv(do.call(rbind,curves),file.path(out,'curves.csv'),row.names=FALSE)
write.csv(do.call(rbind,contrasts),file.path(out,'contrasts.csv'),row.names=FALSE)
capture.output(sessionInfo(),file=file.path(out,'session_info.txt'))
