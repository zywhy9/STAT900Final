## Plot Densities of Prior, Likelihood and Posterior
plotdens <- function(n,y,a,b,beta,title=""){
  curve(dgamma(x,shape=a,rate=b),col="green",xlab=expression(lambda), ylab="Density", main=title, lwd=2, ylim=c(0,1.5), xlim=c(0,5))
  curve(dgamma(x,shape=(a+n),rate=(b+sum(y^beta))),add=T,col="blue", lwd=2)
  # curve(dexp(y,x),add=T,col="red", lwd=2)
  # legend("topright",c("Likelihood","Prior","Posterior"),col=c("red","green","blue"),lwd=rep(2,3),lty=rep(1,3))
  legend("topright",c("Prior","Posterior"),col=c("green","blue"),lwd=rep(2,2),lty=rep(1,2))
}

## Compute the empirical interval length
emlength <- function(x){
  srt <- sort(x)
  num <- length(x)
  lb <- ceiling(num*0.025)
  ub <- ceiling(num*0.975)
  leng <- round(srt[ub]-srt[lb],digits=3)
  return(leng)
}

## Vectorize multiple outcomes
vecx <- function(x,col=1200){
  total <- matrix(NA,nrow=1000,ncol=col)
  for(i in 1:1000){
    Nchain <- nchain(x[[i]])
    xmat <- lapply(x[[i]], as.vector)
    temp <- c()
    for(j in 1:Nchain){
      temp <- c(temp, xmat[[j]])
    }
    total[i,] <- temp
  }
  return(total)
}

## Preplot
"interval.preplot" <- function (x, bin.width = bin.width, max.bins = max.bins,
                                confidence = confidence, autoburnin = autoburnin){
    x <- as.mcmc.list(x)
    nbin <- min(floor((niter(x) - 50)/thin(x)), max.bins)
    if (nbin < 1) {
      stop("Insufficient iterations to produce Interval-ratio plot")
    }
    binw <- floor((niter(x) - 50)/nbin)
    last.iter <- c(seq(from = start(x) + 50 * thin(x), by = binw * 
                         thin(x), length = nbin), end(x))
    shrink <- array(dim = c(nbin + 1, nvar(x), 1))
    dimnames(shrink) <- list(last.iter, varnames(x), "median")
    for (i in 1:(nbin + 1)) {
      shrink[i, , ] <- interval.diag(window(x, end = last.iter[i]), 
                                   confidence = confidence,
                                   autoburnin = autoburnin,
                                   multivariate = FALSE)$psrf[,1]
    }
    all.na <- apply(is.na(shrink[, , 1, drop = FALSE]), 2, all)
    if (any(all.na)) {
      cat("\n******* Error: *******\n")
      cat("Cannot compute Brooks & Gelman's diagnostic for any chain \n")
      cat("segments for variables", varnames(x)[all.na], "\n")
      cat("This indicates convergence failure\n")
    }
    return(list(shrink = shrink, last.iter = last.iter))
  }

## Interval Plot
"interval.plot" <- function (x, bin.width = 10, max.bins = 50, confidence = 0.95,
                            autoburnin = TRUE, auto.layout = TRUE, ask,
                            col = 1:2, lty = 1:2, xlab = "last iteration in chain",
                            ylab = "shrink factor", type = "l", ...){
    if (missing(ask)) {
      ask <- if (is.R()) {
        dev.interactive()
      }
      else {
        interactive()
      }
    }
    x <- as.mcmc.list(x)
    oldpar <- NULL
    on.exit(par(oldpar))
    if (auto.layout) 
      oldpar <- par(mfrow = set.mfrow(Nchains = nchain(x), Nparms = nvar(x)))
    y <- interval.preplot(x, bin.width = bin.width, max.bins = max.bins, 
                        confidence = confidence, autoburnin = autoburnin)
    all.na <- apply(is.na(y$shrink[, , 1, drop = FALSE]), 2, all)
    if (!any(all.na)) 
      for (j in 1:nvar(x)) {
        matplot(y$last.iter, y$shrink[, j, ], col = col, 
                lty = lty, xlab = xlab, ylab = ylab, type = type, 
                ...)
        abline(h = 1)
        ymax <- max(c(1, y$shrink[, j, ]), na.rm = TRUE)
        leg <- dimnames(y$shrink)[[3]]
        # xmax <- max(y$last.iter)
        # legend(xmax, ymax, legend = leg, lty = lty, bty = "n", 
        #        col = col, xjust = 1, yjust = 1)
        title(main = varnames(x)[j])
        if (j==1)
          oldpar <- c(oldpar, par(ask = ask))
      }
    return(invisible(y))
  }

## Interval Diagnostic
"interval.diag" <- function (x, confidence = 0.95, autoburnin=TRUE, multivariate=TRUE){
  x <- as.mcmc.list(x)
  if (nchain(x) < 2) 
    stop("You need at least two chains")
  if (autoburnin && start(x) < end(x)/2 ) 
    x <- window(x, start = end(x)/2 + 1)
  Niter <- niter(x)
  Nchain <- nchain(x)
  Nvar <- nvar(x)
  xnames <- varnames(x)

  
  xmat <- lapply(x, as.matrix)
  LW <- array(sapply(xmat, emlength, simplify=TRUE), dim=c(Nvar,Nvar,Nchain))
  Wmean <- mean(LW)
  total <- c()
  for(i in 1:Nchain){
    total <- c(total, as.vector(xmat[[i]]))
  }
  LT <- emlength(total)
  mpsrf <- NULL
  
  R2.estimate <- LT/Wmean
  R2.upper <- R2.estimate
  psrf <- cbind(R2.estimate, R2.upper)
  dimnames(psrf) <- list(xnames, c("Point est.", "Upper C.I."))
  
  out <- list(psrf = psrf, mpsrf=mpsrf)
  class(out) <- "gelman.diag"
  out
}

"set.mfrow" <-   function (Nchains = 1, Nparms = 1, nplots = 1, sepplot = FALSE) {
    ## Set up dimensions of graphics window: 
    ## If only density plots OR trace plots are requested, dimensions are: 
    ##	1 x 1	if Nparms = 1 
    ##	1 X 2 	if Nparms = 2 
    ##	2 X 2 	if Nparms = 3 or 4 
    ##	3 X 2 	if Nparms = 5 or 6 or 10 - 12 
    ##	3 X 3 	if Nparms = 7 - 9 or >= 13 
    ## If both density plots AND trace plots are requested, dimensions are: 
    ##	1 x 2	if Nparms = 1 
    ##	2 X 2 	if Nparms = 2 
    ##	3 X 2 	if Nparms = 3, 5, 6, 10, 11, or 12 
    ##	4 x 2	if Nparms otherwise 
    ## If separate plots are requested for each chain, dimensions are: 
    ##	1 x 2	if Nparms = 1 & Nchains = 2 
    ##	2 X 2 	if Nparms = 2 & Nchains = 2 OR Nparms = 1 & Nchains = 3 or 4 
    ##	3 x 2	if Nparms = 3 or >= 5 & Nchains = 2  
    ##		   OR Nchains = 5 or 6 or 10 - 12 (and any Nparms) 
    ##	2 x 3	if Nparms = 2 or 4 & Nchains = 3 
    ##	4 x 2   if Nparms = 4 & Nchains = 2  
    ##		   OR Nchains = 4 & Nparms > 1 
    ##	3 x 3	if Nparms = 3 or >= 5  & Nchains = 3  
    ##		   OR Nchains = 7 - 9 or >= 13 (and any Nparms)
    mfrow <- if (sepplot && Nchains > 1 && nplots == 1) {
      ## Separate plots per chain
      ## Only one plot per variable
      if (Nchains == 2) {
        switch(min(Nparms, 5),
               c(1,2),
               c(2,2),
               c(3,2),
               c(4,2),
               c(3,2))
      }
      else if (Nchains == 3) {
        switch(min(Nparms, 5),
               c(2,2),
               c(2,3),
               c(3,3),
               c(2,3),
               c(3,3))
      }
      else if (Nchains == 4) {
        if (Nparms == 1)
          c(2,2)
        else
          c(4,2)
      }
      else if (any(Nchains == c(5,6,10,11,12)))
        c(3,2)
      else if (any(Nchains == c(7,8,9)) || Nchains >=13)
        c(3,3)
      
    }
    else {
      if (nplots==1) {
        ## One plot per variable
        mfrow <- switch(min(Nparms,13),
                        c(1,1),
                        c(1,2),
                        c(2,2),
                        c(2,2),
                        c(3,2),
                        c(3,2),
                        c(3,3),
                        c(3,3),
                        c(3,3),
                        c(3,2),
                        c(3,2),
                        c(3,2),
                        c(3,3))
      }
      else {
        ## Two plot per variable
        ##
        mfrow <- switch(min(Nparms, 13),
                        c(1,2),
                        c(2,2),
                        c(3,2),
                        c(4,2),
                        c(3,2),
                        c(3,2),
                        c(4,2),
                        c(4,2),
                        c(4,2),
                        c(3,2),
                        c(3,2),
                        c(3,2),
                        c(4,2))
      }
    }
    return(mfrow)
}

## JAGS part
jagsres <- function(n,y,a,b,beta,inits,iter=200,thin=1,burnin=0,chain=4){
  start.time <- Sys.time()
  model_string <- "model{
    for(i in 1:n){
      y[i] ~ dweib(beta,lambda) # Model the data
    }
    lambda ~ dgamma(a, b) # The prior
  }"
  jags.data <- list(n=n, y=y, a=a, b=b, beta=beta)
  jags.inits <- list()
  for(i in 1:chain){
    jags.inits[[i]] <- list(lambda=inits[i])
  }
  jags.param <- c("lambda") ## the parameters of the model
  jagsfit <- jags(data=jags.data, n.chains=chain, inits=jags.inits,
                  parameters.to.save=jags.param, n.iter=iter, n.burnin=burnin,n.thin=thin,
                  DIC=FALSE, model.file=textConnection(model_string))
  end.time <- Sys.time()
  time.taken <- as.numeric(end.time - start.time)
  jagsfit.mcmc <- as.mcmc(jagsfit)
  return(list(jagsfit.mcmc,time.taken))
}


## performance
performance <- function(x,a,b,n,y,beta){
  total <- vecx(x)
  posmean <- rep(NA,1000)
  possd <- rep(NA,1000)
  bias <- rep(NA,1000)
  mse <- rep(NA,1000)
  for(i in 1:1000){
    posmean[i] <- unname(summary(x[[i]])$statistics["Mean"])
    possd[i] <- unname(summary(x[[i]])$statistics["SD"])
    
    Nchain <- nchain(x[[i]])
    xmat <- lapply(x[[i]], as.vector)
    total <- c()
    for(j in 1:Nchain){
      total <- c(total, xmat[[j]])
    }
    
    mse[i] <- mean((total-((a+n)/(b+sum(y^beta))))^2)
    bias[i] <- posmean[i]-(a+n)/(b+sum(y^beta))
  }
  
  res <- list(posteriormean=mean(posmean),bias=mean(bias),mse=mean(mse),posteriorsd=mean(possd))
  return(res)
}


## Within Chain MCSE
within_chain_se <- function(x,realv,type="bias"){
  ssw <- 0
  if(type=="bias"){
    for(i in 1:length(x)){
      ssw <- ssw + sum((x[[i]]-mean(x[[i]]))^2)
    }
  }else if(type=="mse"){
    for(i in 1:length(x)){
      msei <- mean((x[[i]]-realv)^2)
      uniti <- (x[[i]]-realv)^2
      ssw <- ssw + sum((uniti-msei)^2)
    }
  }
  ssw <- sqrt(ssw/(length(x)*length(x[[1]])*(length(x)*length(x[[1]])-1)))
  return(ssw)
}

## Between Chain MCSE
between_chain_se <- function(x,realv,type="bias"){
  ssb <- 0
  ntem <- length(x[[1]])
  if(type=="bias"){
    meant <- unname(summary(x)$statistics["Mean"])
    for(i in 1:length(x)){
      meani <- mean(x[[i]])
      ssb <- ssb + ntem * (meani-meant)^2
    }
  }else if(type=="mse"){
    Nchain <- nchain(x)
    xmat <- lapply(x, as.vector)
    total <- c()
    for(j in 1:Nchain){
      total <- c(total, xmat[[j]])
    }
    mset <- mean((total-realv)^2)
    for(i in 1:length(x)){
      msei <- mean((x[[i]]-realv)^2)
      ssb <- ssb + ntem*(msei-mset)^2
    }
  }
  ssb <- sqrt(ssb/(length(x)*ntem*(length(x)*ntem-1)))
  return(ssb)
}

## Total Chain MCSE
total_chain_se <- function(x,realv,type="bias"){
  sst <- 0
  if(type=="bias"){
    meant <- unname(summary(x)$statistics["Mean"])
    for(i in 1:length(x)){
      sst <- sst + sum((x[[i]]-meant)^2)
    }
  }else if(type=="mse"){
    Nchain <- nchain(x)
    xmat <- lapply(x, as.vector)
    total <- c()
    for(j in 1:Nchain){
      total <- c(total, xmat[[j]])
    }
    mset <- mean((total-realv)^2)
    for(i in 1:length(x)){
      uniti <- (x[[i]]-realv)^2
      sst <- sst + sum((uniti-mset)^2)
    }
  }
  sst <- sqrt(sst/(length(x)*length(x[[1]])*(length(x)*length(x[[1]])-1)))
  return(sst)
}

## Within Simulation MCSE
within_sim_se <- function(x,realv,type="bias"){
  ssw <- 0
  nt <- length(x[[1]])*length(x[[1]][[1]])
  total <- vecx(x,col=nt)
  if(type=="bias"){
    for(i in 1:1000){
      ssw <- ssw + sum((total[i,]-mean(total[i,]))^2)
    }
  }else if(type=="mse"){
    for(i in 1:1000){
      msei <- mean((total[i,]-realv)^2)
      uniti <- (total[i,]-realv)^2
      ssw <- ssw + sum((uniti-msei)^2)
    }
  }
  ssw <- sqrt(ssw/(nt*1000*(nt*1000-1)))
  return(ssw)
}

## Between Simulation MCSE
between_sim_se <- function(x,realv,type="bias"){
  ssb <- 0
  nt <- length(x[[1]])*length(x[[1]][[1]])
  total <- vecx(x,col=nt)
  if(type=="bias"){
    meant <- mean(total)
    for(i in 1:1000){
      meani <- mean(total[i,])
      ssb <- ssb + nt * (meani-meant)^2
    }
  }else if(type=="mse"){
    mset <- mean((total-realv)^2)
    for(i in 1:1000){
      msei <- mean((total[i,]-realv)^2)
      ssb <- ssb + nt*(msei-mset)^2
    }
  }
  ssb <- sqrt(ssb/(nt*1000*(nt*1000-1)))
  return(ssb)
}

## Total Simulation MCSE
total_sim_se <- function(x,realv,type="bias"){
  sst <- 0
  nt <- length(x[[1]])*length(x[[1]][[1]])
  total <- vecx(x,col=nt)
  if(type=="bias"){
    meant <- mean(total)
    for(i in 1:1000){
      sst <- sst + sum((total[i,]-meant)^2)
    }
  }else if(type=="mse"){
    mset <- mean((total-realv)^2)
    for(i in 1:1000){
      uniti <- (total[i,]-realv)^2
      sst <- sst + sum((uniti-mset)^2)
    }
  }
  sst <- sqrt(sst/(nt*1000*(nt*1000-1)))
  return(sst)
}
