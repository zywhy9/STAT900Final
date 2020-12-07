library(R2jags)
library(coda)
library(MCMCpack)
source("function.R")

set.seed(1)

## Set Model
n <- 10
beta <- 1.5
theta <- 1
y <- t(matrix(rweibull(n*10,beta,theta),ncol=n))
# curve(dweibull(x,beta,theta),col="red", lwd=2, xlim=c(0,5), main="Likelihood Function")
a <- 2
b <- 2

## Model
saveadd <- "simulation/6x190/"

unittime <- matrix(0,nrow=10,ncol=1000)
colnames(unittime) <- 1:1000

start.time0 <- Sys.time()
time1 <- rep(0,10)
for(i in 1:10){
  start.time1 <- Sys.time()
  for(j in 1:1000){
    out <- jagsres(n,y[i,],a,b,beta,inits=rgamma(6,shape=(a+n),rate=(b+sum(y^beta))),iter=190,chain=6)
    unittime[i,j] <- out[[2]]
    saveRDS(out,paste0(saveadd,"out",i,"-",j))
    rm(out)
  }
  end.time1 <- Sys.time()
  time1[i] <- end.time1 - start.time1
}
end.time0 <- Sys.time()
time0 <- end.time0 - start.time0

write.csv(unittime,"simulation/6x190time.csv")