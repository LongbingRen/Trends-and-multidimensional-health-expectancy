fit_dimension <- function(d,dimension,cfg,extra_covariates=character()) {
  z <- d; z$state <- z[[dimension]]
  Q <- rbind(c(0,.04,.03),c(.08,0,.08),c(0,0,0))
  if(dimension=='disease') Q[2,1] <- 0
  # Living observations are panel censored. Only verified exact deaths are 3.
  z$obs_type <- ifelse(z$dead==1 & z$death_exact==1,3L,1L)
  covars <- c(cfg$model_covariates,extra_covariates)
  if(anyNA(z[c('state','age',covars)])) stop('Missing model variables.')
  design <- model.matrix(reformulate(covars),z)
  if(qr(design)$rank<ncol(design)) stop('Rank-deficient covariate design: ',dimension)
  fit <- msm::msm(state~age,subject=sid,data=z,qmatrix=Q,
    covariates=reformulate(covars),center=FALSE,obstype=obs_type,
    method='BFGS',control=list(maxit=10000,reltol=1e-9))
  if(!is.null(fit$opt$convergence) && fit$opt$convergence!=0)
    stop('Model failed to converge: ',dimension)
  if(is.null(fit$covmat) || any(!is.finite(fit$covmat)))
    stop('Invalid parameter covariance: ',dimension)
  if(min(eigen(fit$covmat,symmetric=TRUE,only.values=TRUE)$values)<=0)
    stop('Non-positive-definite parameter covariance: ',dimension)
  fit
}
fit_all <- function(d,cfg) {
  out <- lapply(cfg$dimensions,function(dim) {
    message('Fitting ',unique(d$age_group),': ',dim)
    fit_dimension(d,dim,cfg)
  })
  names(out) <- cfg$dimensions; out
}
baseline_models <- function(d,cfg) {
  b <- d[d$baseline==1,]
  followup <- do.call(rbind,lapply(split(d,d$sid),function(z) {
    data.frame(sid=z$sid[1],time=as.numeric(tail(z$date,1)-z$date[1])/365.25,
      event=tail(z$dead,1))
  }))
  b <- merge(b,followup,by='sid',sort=FALSE)
  demography <- c('period','age_actual','female','urban','educated','with_others','married')
  mortality <- c(demography,'n_conditions','physical','cognition','care','srh')
  rows <- list(); objects <- list()
  for(a in cfg$centres) {
    z <- b[b$age_group==a,]
    if(anyNA(z[c(mortality,'time','event')])) stop('Missing baseline model data.')
    cox <- survival::coxph(reformulate(mortality,'survival::Surv(time,event)'),
      data=z,ties='efron',x=TRUE)
    objects[[paste0('cox_',a)]] <- cox
    sc <- summary(cox)$coefficients['period',]
    rows[[length(rows)+1L]] <- data.frame(age=a,outcome='mortality',type='HR',
      estimate=exp(sc['coef']),lower=exp(sc['coef']-1.96*sc['se(coef)']),
      upper=exp(sc['coef']+1.96*sc['se(coef)']),p=sc['Pr(>|z|)'])
    objects[[paste0('km_',a)]] <- survival::survfit(survival::Surv(time,event)~period,data=z)
    objects[[paste0('logrank_',a)]] <- survival::survdiff(survival::Surv(time,event)~period,data=z)
    objects[[paste0('ph_',a)]] <- survival::cox.zph(cox)
    for(dim in cfg$dimensions) {
      z$adverse <- as.integer(z[[dim]]==2)
      fit <- glm(reformulate(demography,'adverse'),data=z,family=binomial())
      if(!fit$converged) stop('Logistic model failed: ',dim)
      s <- summary(fit)$coefficients['period',]
      rows[[length(rows)+1L]] <- data.frame(age=a,outcome=dim,type='OR',
        estimate=exp(s[1]),lower=exp(s[1]-1.96*s[2]),
        upper=exp(s[1]+1.96*s[2]),p=s[4])
      objects[[paste(dim,a,sep='_')]] <- fit
    }
  }
  list(estimates=do.call(rbind,rows),models=objects)
}
distribution_table <- function(d,cfg,weighted=FALSE) {
  b <- d[d$baseline==1,]
  vars <- cfg$decomposition_factors
  out <- list()
  for(a in cfg$centres) for(p in 0:1) {
    z <- b[b$age_group==a & b$period==p,]
    w <- if(weighted) z$survey_weight else rep(1,nrow(z))
    if(anyNA(w) || any(w<0) || sum(w)<=0) stop('Invalid baseline weights.')
    for(v in vars) out[[length(out)+1L]] <- data.frame(age=a,period=p,
      factor=v,prevalence=weighted.mean(z[[v]],w),n=nrow(z),weighted=weighted)
  }
  do.call(rbind,out)
}
