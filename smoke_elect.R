# Interface smoke test using the ELECT package's own simulated example data.
# Run after installation: Rscript tests/smoke_elect.R
source('config.R'); source('R/03_life_expectancy.R')
data('electData',package='elect')
Q <- rbind(c(0,.01,.01),c(0,0,.01),c(0,0,0))
fit <- msm::msm(state~age,subject=id,data=electData,qmatrix=Q,
  center=FALSE,deathexact=3,covariates=~age+x)
cfg$model_covariates <- c('age','x'); cfg$max_age <- 120
cfg$S <- 4L; cfg$h <- .5
base <- electData[electData$bsline==1,]
a <- elect_profile(fit,list(age=0,x=1),base,70,cfg,123)
b <- elect_profile(fit,list(age=0,x=1),base,70,cfg,123)
stopifnot(all(dim(a)==c(5,3)),isTRUE(all.equal(a,b)),all(is.finite(a)))
cat('ELECT interface smoke test passed. Not a CLHLS reproduction test.\n')
