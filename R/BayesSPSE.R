BayesSPSE <- function(x,y,nk){
  require(jagsUI)
  
  ### Assemble data into list for JAGS====
  y     = as.vector(scale(y))
  x     = x
  n     = length(x)
  nk    = tryCatch(nk, error=function(e) (round(sqrt(n-2),0) + 1) )
  knots = as.vector(c(0,quantile(x, seq(0,1,len=nk)),1))
  dataList = list( y=y, x=x, n=n, nk=nk, knots=knots )
  
  # Define the model====
  modelString = "
  model {
    ### Frontier estimation
    for ( i in 1:n ) {
      y[i]    ~ dnorm( pred[i] , err )
      pred[i] <- b0 + sum(Beta[1:nk] * Z[i,1:nk]) + v[i]
      v[i]    ~ dexp( parF )
      eff[i]  <- 1/exp(v[i])
    }
    err  ~ dgamma( shape , rate )
    
    ### Additive model
    b0 ~ dnorm( mub0 , sgm0 )
    for (l in 1:nk) {
      Beta[l] ~ ddexp( mubt , tau )
      dgr[l]  <- dg[l] + 1
      dg[l]   ~ dbinom( theta , l )
      for (i in 1:n) {
        Z[i,l] <- ((x[i]^dgr[l]) - knots[l]) / sd(x)
      }
    }

    ### Overall priors
    # Error priors
    shape  ~ dunif( 1e-2, 1e2 )
    rate   ~ dunif( 1e-2, 1e2 )

    # Efficiency priors
    parF  ~ dgamma( 1e-2, 1e-2 )

    # Regression priors
    mubt   ~ dnorm( 0 , 1 )
    tau    ~ dunif( 1e-2, 1e2 )
    mub0   ~ dnorm( 0 , 1 )
    sgm0   ~ dunif( 1e-2, 1e2 )
    theta  ~ dbeta( 1 , 1 )
    
  }
  " # close quote for modelString
  model = textConnection(modelString)
  
  # Run the chains====
  # Name the parameters to be monitored
  params <- c("eff","pred","Beta",'b0','dgr')
  # Define some MCMC parameters for JAGS
  nthin    = 1    # How Much Thinning?
  nchains  = 4    # How Many Chains?
  nburnin  = 200  # How Many Burn-in Samples?
  nadapt   = 500  # How Many adaptation Samples?
  nsamples = 1700 # How Many Recorded Samples?
  ### Calling JAGS to sample
  startTime = proc.time()
  set.seed(666)
  samples <- jagsUI::jags(dataList, NULL, params, model.file=model,
                          n.chains=nchains, n.adapt=nadapt, n.iter=nsamples, 
                          n.burnin=nburnin, n.thin=nthin, DIC=T)
  stopTime = proc.time(); elapsedTime = stopTime - startTime; methods::show(elapsedTime)
  
  ### Inspect and diagnose the run====
  GRD <- tryCatch(coda::gelman.diag(samples$samples), error=function(e) NA)
  eff <- colMeans(samples$sims.list$eff)
  pred <- colMeans(samples$sims.list$pred)
  betas <- colMeans(samples$sims.list$Beta)
  b0    <- mean(samples$sims.list$b0)
  dic <- samples$DIC
  full <- samples
  Result <- list("eff"=eff,"pred"=pred,"GRD"=GRD,'betas'=betas,
                 "b0"=b0,"dic"=dic,"full"=full)
  return(Result)
}
