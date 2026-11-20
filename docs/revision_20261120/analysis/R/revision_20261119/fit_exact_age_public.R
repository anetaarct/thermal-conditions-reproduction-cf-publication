.libPaths(c(path.expand('~/R/cf_model_comparison_20261013_R46'),
            path.expand('~/R/x86_64-pc-linux-gnu-library/4.5'), .libPaths()))
Sys.setenv(OMP_NUM_THREADS = 1, OPENBLAS_NUM_THREADS = 1)
library(glmmTMB)
out <- 'age_sensitivity_20261119'
dir.create(out, showWarnings = FALSE)
paths <- c(clutch = 'data_derived/analysis_clean_20261119/clutch_model_data.csv',
           fledging = 'data_derived/analysis_clean_20261119/fledging_model_data.csv')
models <- c(clutch='models/clean_20261119/clutch_baseline_selected.rds', fledging='models/clean_20261119/fledging_baseline_selected.rds')
fixed <- diagnostics <- samples <- age_counts <- list()
for (trait in names(paths)) {
  d <- read.csv(paths[[trait]])
  known <- read.csv(file.path(out,paste0(trait,'_known_age_model_data.csv')))
  keys <- paste(d$female_id,d$year,sep='|')
  known_keys <- paste(known$female_id,known$year,sep='|')
  stopifnot(!anyDuplicated(keys),!anyDuplicated(known_keys),all(known_keys %in% keys))
  d$exact_age <- known$exact_age[match(keys,known_keys)]
  d$female_id <- factor(d$female_id)
  d$year_f <- factor(d$year)
  d$age_class <- factor(d$age_class,levels=c('YOUNG','OLD'))
  d$temp_between <- ave(d$temp_mean_c, d$female_id, FUN = mean)
  d$temp_within <- d$temp_mean_c - d$temp_between
  known <- !is.na(d$exact_age)
  cat(trait, 'full', nrow(d), 'known', sum(known), 'females',
      length(unique(d$female_id[known])), '\n'); flush.console()
  samples[[trait]] <- data.frame(trait, records_full = nrow(d), females_full = length(unique(d$female_id)),
    records_known = sum(known), females_known = length(unique(d$female_id[known])), seasons_known = length(unique(d$year_f[known])))
  age_counts[[trait]] <- data.frame(trait, age = names(table(d$exact_age)), records = as.integer(table(d$exact_age)))
  write.csv(do.call(rbind, samples), file.path(out, 'sample_counts.csv'), row.names = FALSE)
  write.csv(do.call(rbind, age_counts), file.path(out, 'age_counts.csv'), row.names = FALSE)
  d <- droplevels(d[known, ])
  stopifnot(all(d$exact_age >= 1), all(d$exact_age[d$age_class == 'YOUNG'] == 1), all(d$exact_age[d$age_class == 'OLD'] >= 2))
  d$age_c <- d$exact_age - mean(d$exact_age)
  d$age_c2 <- d$age_c^2
  old <- readRDS(models[[trait]])
  for (kind in c('age_class', 'exact_age_linear', 'exact_age_quadratic')) {
    age_terms <- switch(kind, age_class = 'age_class', exact_age_linear = 'age_c', exact_age_quadratic = 'age_c + age_c2')
    rhs <- paste('rel_LD + temp_within + temp_between + year_c +', age_terms)
    response <- if (trait == 'clutch') 'clutch' else 'cbind(fledged, clutch - fledged)'
    formula <- as.formula(paste(response, '~', rhs, '+ (1|female_id) + (1|year_f)'))
    disp <- as.formula(paste('~', rhs, if (trait == 'fledging') '+ (1|year_f)' else ''))
    reference_file <- file.path(out,paste0(trait,'_age_class_refined.rds'))
    reference <- if(file.exists(reference_file))readRDS(reference_file) else old
    class_file<-file.path(out,paste0(trait,'_age_class.rds'))
    if(!file.exists(reference_file)&&file.exists(class_file)) {
      class_fit<-readRDS(class_file)
      if(isTRUE(class_fit$sdr$pdHess)&&class_fit$fit$convergence==0)reference<-class_fit
    }
    if(kind=='exact_age_quadratic') {
      linear_file<-file.path(out,paste0(trait,'_exact_age_linear.rds'))
      if(file.exists(linear_file)) {
        candidate_reference<-readRDS(linear_file)
        if(isTRUE(candidate_reference$sdr$pdHess)&&candidate_reference$fit$convergence==0)reference<-candidate_reference
      }
    }
    cf <- fixef(reference)
    make_start <- function(part) {
      cols <- colnames(model.matrix(as.formula(paste('~', rhs)), d))
      b <- setNames(rep(0, length(cols)), cols)
      common <- intersect(cols, names(cf[[part]]))
      b[common] <- cf[[part]][common]
      if('temp_mean_c' %in% names(cf[[part]]))b[c('temp_within', 'temp_between')] <- cf[[part]]['temp_mean_c']
      if(kind!='age_class' && 'age_classOLD' %in% names(cf[[part]])) {
        age_regression<-if(kind=='exact_age_quadratic')lm(as.numeric(age_class=='OLD')~age_c+age_c2,data=d) else lm(as.numeric(age_class=='OLD')~age_c,data=d)
        projection<-coef(age_regression)
        b['(Intercept)']<-b['(Intercept)']+cf[[part]]['age_classOLD']*projection['(Intercept)']
        b['age_c']<-cf[[part]]['age_classOLD']*projection['age_c']
        if(kind=='exact_age_quadratic')b['age_c2']<-cf[[part]]['age_classOLD']*projection['age_c2']
      }
      unname(b)
    }
    start <- list(beta = make_start('cond'), betadisp = make_start('disp'))
    if (trait == 'fledging') start$betazi <- unname(cf$zi)
    parameters<-reference$obj$env$parList(par=reference$fit$parfull)
    for(name in intersect(c('theta','thetadisp','thetazi'),names(parameters)))if(length(parameters[[name]]))start[[name]]<-parameters[[name]]
    age_optimizer<-function(par,fn,gr,control=list(),...) {
      zi_names<-if(length(cf$zi))paste0('zi~',names(cf$zi)) else character(0)
      old_names<-c(paste0('cond~',names(cf$cond)),zi_names,paste0('disp~',names(cf$disp)))
      old_names<-sub('age_classOLD','age_c',old_names,fixed=TRUE)
      random_count<-nrow(reference$sdr$cov.fixed)-length(old_names)
      old_names<-c(old_names,paste0('random~',seq_len(random_count)))
      columns<-colnames(model.matrix(as.formula(paste('~',rhs)),d))
      new_names<-c(paste0('cond~',columns),zi_names,paste0('disp~',columns),paste0('random~',seq_len(random_count)))
      stopifnot(length(new_names)==length(par))
      indices<-match(new_names,old_names)
      A<-matrix(0,length(new_names),length(old_names))
      matched<-which(!is.na(indices));A[cbind(matched,indices[matched])]<-1
      V<-A%*%reference$sdr$cov.fixed%*%t(A)
      missing<-which(is.na(indices))
      for(j in missing)V[j,j]<-if(grepl('^cond',new_names[j])) .005^2 else .05^2
      eig<-eigen((V+t(V))/2,symmetric=TRUE);stopifnot(min(eig$values)>0)
      S<-eig$vectors%*%diag(sqrt(eig$values))
      if(trait=='clutch') {
        inner<-environment(fn)$inner.control;inner$maxit<-100L;inner$grad.tol<-1e-7
        environment(fn)$inner.control<-inner
      }
      ff<-function(q)fn(drop(par+S%*%q));gg<-function(q)drop(crossprod(S,gr(drop(par+S%*%q))))
      z<-optim(rep(0,length(par)),ff,gg,method='BFGS',control=list(maxit=1500,reltol=1e-12,trace=1,REPORT=5))
      list(par=drop(par+S%*%z$par),objective=z$value,convergence=z$convergence,message=z$message)
    }
    file <- file.path(out, paste0(trait, '_', kind, '.rds'))
    cat(format(Sys.time()), trait, kind, 'START\n'); flush.console()
    if(file.exists(file)) {
      cached<-readRDS(file)
      if(!isTRUE(cached$sdr$pdHess) || cached$fit$convergence!=0) {
        archive<-sub('\\.rds$',paste0('_inadmissible_',format(Sys.time(),'%H%M%S'),'.rds'),file)
        stopifnot(file.rename(file,archive))
      }
    }
    if (file.exists(file)) m <- readRDS(file) else {
      m <- glmmTMB(formula, dispformula = disp,
        ziformula = if (trait == 'clutch') ~0 else ~1 + (1|year_f),
        family = if (trait == 'clutch') compois() else betabinomial(), data = d,
        start = start, control = if(kind=='age_class')glmmTMBControl(parallel = 4,
          optCtrl = list(iter.max = 3000, eval.max = 6000)) else glmmTMBControl(optimizer=age_optimizer,parallel=4))
      saveRDS(m, file)
    }
    dg <- data.frame(trait, model = kind, n = nobs(m), convergence = m$fit$convergence,
      pdHess = isTRUE(m$sdr$pdHess), logLik = as.numeric(logLik(m)), AIC = AIC(m),
      scaled_gradient = max(abs(m$sdr$gradient.fixed * sqrt(diag(m$sdr$cov.fixed)))))
    diagnostics[[paste(trait, kind)]] <- dg
    write.csv(do.call(rbind, diagnostics), file.path(out, 'diagnostics.csv'), row.names = FALSE)
    print(dg); flush.console()
    if (dg$convergence == 0 && dg$pdHess && is.finite(dg$scaled_gradient) && dg$scaled_gradient < .005) {
      co <- summary(m)$coefficients
      fixed[[paste(trait, kind)]] <- do.call(rbind, lapply(names(co), function(component) {
        x <- co[[component]]
        if (!length(x)) return(NULL)
        data.frame(trait, model = kind, component, term = rownames(x), estimate = x[,1], SE = x[,2],
          lower = x[,1] - qnorm(.975)*x[,2], upper = x[,1] + qnorm(.975)*x[,2])
      }))
      write.csv(do.call(rbind, fixed), file.path(out, 'fixed_effects.csv'), row.names = FALSE)
    }
    capture.output(summary(m), file = file.path(out, paste0(trait, '_', kind, '_summary.txt')))
  }
}
capture.output(sessionInfo(), file = file.path(out, 'session_info.txt'))
cat('COMPLETE\n')
