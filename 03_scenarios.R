run_stratum <- function(d,cfg,outdir) {
  centre <- unique(d$age_group); stopifnot(length(centre)==1L)
  models <- fit_all(d,cfg)
  saveRDS(models,file.path(outdir,paste0('models_',centre,'.rds')))
  hr <- lapply(models,msm::hazard.msm)
  capture.output(hr,file=file.path(outdir,paste0('transition_HR_',centre,'.txt')))
  p <- lapply(0:1,function(k) profiles(d[d$period==k,],cfg))
  seed <- cfg$seed+as.integer(centre)*10000L
  base <- lapply(0:1,function(k) population_le(models,d,p[[k+1]],k,centre,cfg,
                                             seed+k*100L))
  # Preserve all underlying scenarios as well as final contrasts.
  scenarios <- list(period0=base[[1]],period1=base[[2]])
  rows <- list()
  add <- function(x,label,period=NA_integer_,factor=NA_character_,
                  target=NA_real_,sex='all') {
    rows[[length(rows)+1L]] <<- summarise_le(x,centre,label,period,factor,target,sex)
  }
  for(k in 0:1) {
    add(base[[k+1]]$raw,'raw',k)
    add(base[[k+1]]$harmonised,'harmonised',k)
    add(healthy_proportion(base[[k+1]]$harmonised),'proportion',k)
  }
  difference <- base[[2]]$harmonised-base[[1]]$harmonised
  add(difference,'period_difference')
  if(cfg$run_sex_specific) for(s in 0:1) {
    sx <- vector('list',2)
    for(k in 0:1) {
      ps <- p[[k+1]][p[[k+1]]$female==s,]
      if(!nrow(ps)) stop('No participants for sex-specific estimate.')
      ps$weight <- ps$weight/sum(ps$weight)
      sx[[k+1]] <- population_le(models,d,ps,k,centre,cfg,
        seed+2000L+s*100L+k*10L,sex=s)
      add(sx[[k+1]]$harmonised,'harmonised',k,sex=ifelse(s==1,'women','men'))
    }
    add(sx[[2]]$harmonised-sx[[1]]$harmonised,'period_difference',
      sex=ifelse(s==1,'women','men'))
    scenarios[[paste0('sex_',s)]] <- sx
  }
  if(cfg$run_decomposition) for(j in seq_along(cfg$decomposition_factors)) {
    f <- cfg$decomposition_factors[j]
    message('Decomposition ',centre,': ',f)
    p21 <- transport_factor(p[[2]],p[[1]],f)
    p12 <- transport_factor(p[[1]],p[[2]],f)
    L21 <- population_le(models,d,p21,1,centre,cfg,seed+10000L+j*100L)
    L12 <- population_le(models,d,p12,0,centre,cfg,seed+10000L+j*100L+10L)
    contribution <- symmetric_contribution(base[[1]]$harmonised,
      base[[2]]$harmonised,L21$harmonised,L12$harmonised)
    add(contribution,'decomposition',factor=f)
    scenarios[[paste0('decomposition_',f)]] <- list(L11=base[[1]],
      L22=base[[2]],L21=L21,L12=L12,contribution=contribution)
  }
  if(cfg$run_counterfactual) for(target in c(1,.8)) {
    for(j in seq_along(cfg$counterfactual_factors)) {
      f <- cfg$counterfactual_factors[j]
      pc <- attainment_profiles(p[[2]],f,target)
      if(is.null(pc)) next # Coverage already at / above target: do not reduce it.
      message('Counterfactual ',centre,': ',f,' target ',target)
      cf <- population_le(models,d,pc,1,centre,cfg,
                          seed+20000L+j*100L+as.integer(target*10))
      gain <- cf$harmonised-base[[2]]$harmonised
      add(cf$harmonised,'counterfactual_absolute',1,f,target)
      add(gain,'counterfactual_gain',1,f,target)
      scenarios[[paste0('cf_',f,'_',target)]] <- list(reference=base[[2]],
        counterfactual=cf,gain=gain)
    }
  }
  saveRDS(scenarios,file.path(outdir,paste0('scenario_draws_',centre,'.rds')))
  tab <- do.call(rbind,rows)
  write.csv(tab,file.path(outdir,paste0('estimates_',centre,'.csv')),row.names=FALSE)
  tab
}
