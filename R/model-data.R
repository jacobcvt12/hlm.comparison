DP <- function() {
  # data
  data(indon)
  attach(indon)

  baseage2<-baseage**2
  follow<-age-baseage
  follow2<-follow**2

  # prior information
  beta0 <- rep(0, 9)
  Sbeta0 <- diag(1000, 9)
  tinv <- diag(1, 1)
  prior <- list(a0=2, b0=0.1, nu0=4, tinv=tinv, mub=rep(0, 1),
                Sb=diag(1000, 1), beta0=beta0, Sbeta0=Sbeta0)

  # initial state
  state <- NULL

  # MCMC paramaters
  nburn <- 5000
  nsave <- 5000
  nskip <- 0
  ndisplay <- 1000
  mcmc <- list(nburn=nburn, nsave=nsave, nskip=nskip, ndisplay=ndisplay)

  # fit the model
  fit1 <- DPglmm(fixed=infect~gender+height+cosv+sinv+
                       xero+baseage+baseage2+follow+follow2,
                 random=~1|id,
                 family=binomial(logit),
                 prior=prior,
                 mcmc=mcmc,
                 state=state,
                 status=TRUE)

  return(fit1)
}
