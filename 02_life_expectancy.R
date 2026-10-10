elect_profile <- function(model,covariates,sddata,centre,cfg,seed,S=cfg$S) {
  # The first covariate MUST be age. All covariates numeric, center=FALSE.
  # Both target age and upper age limit use the centred time scale.
  stopifnot(identical(names(covariates),cfg$model_covariates))
  if(!setequal(unique(sddata$state),1:2))
    stop('Both living states must occur in state-distribution data.')
  el <- elect::elect(x=model,b.covariates=covariates,
    statedistdata=sddata,time.scale.msm='years',
    statedist.covariates='age',h=cfg$h,age.max=cfg$max_age-centre,
    S=S,setseed=seed,method=cfg$method)
  idx <- match(c('e.1','e.2','e'),names(el$pnt))
  if(anyNA(idx)) stop('Unexpected ELECT output names; inspect package version.')
  point <- as.numeric(el$pnt[idx])
  sims <- if(S>0L) as.matrix(el$sim)[,idx,drop=FALSE] else matrix(numeric(),0,3)
  result <- rbind(point,sims)
  colnames(result) <- c('healthy','unhealthy','total')
  if(any(!is.finite(result))) stop('Non-finite ELECT results; inspect fit / draws.')
  if(max(abs(result[,1]+result[,2]-result[,3]))>1e-7)
    stop('ELECT healthy + unhealthy != total.')
  if(any(result < -1e-8) || any(result[,3] > cfg$max_age-centre+cfg$h))
    stop('Life expectancy outside the integration bounds.')
  result
}
population_le <- function(models,d,p,period,centre,cfg,seed,sex=NULL) {
  # Cohort-period transition regime is selected through period; no refitting
  # separate period models, consistent with the main text's pooled-period fit.
  vars <- setdiff(cfg$model_covariates,'age')
  p$period <- period
  if(cfg$standardisation=='mean_profile') {
    q <- as.data.frame(lapply(p[vars],weighted.mean,w=p$weight))
    q$weight <- 1
  } else if(cfg$standardisation=='empirical_profiles') {
    q <- compress_profiles(p,vars)
  } else stop('Unknown standardisation.')
  base <- d[d$baseline==1 & d$period==period,]
  if(!is.null(sex)) base <- base[base$female==sex,]
  ans <- array(NA_real_,dim=c(cfg$S+1L,5L,3L),
    dimnames=list(NULL,cfg$dimensions,c('healthy','unhealthy','total')))
  for(j in seq_along(cfg$dimensions)) {
    dim <- cfg$dimensions[j]; sddata <- base; sddata$state <- base[[dim]]
    acc <- matrix(0,cfg$S+1L,3L)
    for(i in seq_len(nrow(q))) {
      covs <- as.list(c(age=0,unlist(q[i,vars,drop=FALSE])))
      covs <- covs[cfg$model_covariates]
      # Same seed across profiles WITHIN one dimension/scenario aligns draws
      # of the same fitted models in ELECT 1.2/1.2.1. Verify after any upgrade.
      e <- elect_profile(models[[dim]],covs,sddata,centre,cfg,seed*10L+j)
      acc <- acc+q$weight[i]*e
    }
    ans[,j,] <- acc
  }
  list(raw=ans,harmonised=harmonise_le(ans))
}
harmonise_le <- function(x) {
  # Apply to each simulation row BEFORE SDs, CIs or scenario contrasts.
  stopifnot(length(dim(x))==3L,dim(x)[2]==5L,dim(x)[3]==3L)
  out <- x
  for(b in seq_len(dim(x)[1])) {
    total <- x[b,,3]
    if(any(total<=0)) stop('Non-positive LE in harmonisation.')
    common <- mean(total)
    out[b,,1] <- x[b,,1]*common/total
    out[b,,2] <- x[b,,2]*common/total
    out[b,,3] <- common
  }
  out
}
summarise_le <- function(x,age,scenario,period=NA_integer_,factor=NA_character_,
                         target=NA_real_,sex='all') {
  out <- list()
  for(j in seq_len(dim(x)[2])) for(k in seq_len(dim(x)[3])) {
    point <- x[1,j,k]
    sims <- if(dim(x)[1]>1) x[-1,j,k] else numeric()
    se <- if(length(sims)>1) sd(sims) else NA_real_
    quant <- if(length(sims)>1) quantile(sims,c(.025,.975),names=FALSE) else c(NA,NA)
    out[[length(out)+1L]] <- data.frame(age=age,scenario=scenario,period=period,
      factor=factor,target=target,sex=sex,dimension=dimnames(x)[[2]][j],
      metric=dimnames(x)[[3]][k],estimate=point,se=se,
      lower=point-1.96*se,upper=point+1.96*se,
      quantile_lower=quant[1],quantile_upper=quant[2])
  }
  do.call(rbind,out)
}
healthy_proportion <- function(x) {
  out <- x[,,1,drop=FALSE]/x[,,3,drop=FALSE]
  dimnames(out)[[3]] <- 'healthy_proportion'; out
}
symmetric_contribution <- function(L11,L22,L21,L12) {
  # L21: period 2 regime / other covariates, factor from period 1.
  # L12: period 1 regime / other covariates, factor from period 2.
  .5*((L12-L11)+(L22-L21))
}
