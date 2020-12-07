library(R2jags)
library(coda)
library(MCMCpack)
source("function.R")

set.seed(1)

output <- matrix(NA,nrow=10,ncol=20)

for(k in 1:10){
  dsno <- k
  ## Read Data
  out <- list()
  for(i in 1:1000){
    out[[i]] <- readRDS(paste0("simulation/6x200/out",dsno,"-",i))[[1]]
  }
  
  ## Summary
  perm <- performance(out,a,b,n,y[dsno,],beta)
  
  # (a+n)/(b+sum(y[dsno,]^beta))
  # sqrt((a+n)/(b+sum(y[dsno,]^beta))^2)
  
  ## MCSE
  res <- matrix(NA,nrow=6,ncol=1000)
  for(i in 1:1000){
    res[1,i] <- within_chain_se(out[[i]],realmean[dsno],type="bias")
    res[2,i] <- between_chain_se(out[[i]],realmean[dsno],type="bias")
    res[3,i] <- total_chain_se(out[[i]],realmean[dsno],type="bias")
    res[4,i] <- within_chain_se(out[[i]],realmean[dsno],type="mse")
    res[5,i] <- between_chain_se(out[[i]],realmean[dsno],type="mse")
    res[6,i] <- total_chain_se(out[[i]],realmean[dsno],type="mse")
  }
  cv6 <- mean(res[1,])
  cv7 <- mean(res[2,])
  cv8 <- mean(res[3,])
  cv13 <- mean(res[4,])
  cv14 <- mean(res[5,])
  cv15 <- mean(res[6,])
  
  cv9 <- within_sim_se(out,realv=realmean[dsno],type="bias")
  cv10 <- between_sim_se(out,realv=realmean[dsno],type="bias")
  cv11 <- total_sim_se(out,realv=realmean[dsno],type="bias")
  cv16 <- within_sim_se(out,realv=realmean[dsno],type="mse")
  cv17 <- between_sim_se(out,realv=realmean[dsno],type="mse")
  cv18 <- total_sim_se(out,realv=realmean[dsno],type="mse")
  
  ## Convergence Diagnostic
  diag <- matrix(NA,nrow=2,ncol=1000)
  for(i in 1:1000){
    diag[1,i] <- gelman.diag(out[[i]])$psrf[1,1]
    diag[2,i] <- interval.diag(out[[i]])$psrf[1,1]
  }
  cv19 <- mean(diag[1,])
  cv20 <- mean(diag[2,])
  
  ## Output
  cv1 <- realmean[dsno]
  cv2 <- unname(perm[1][[1]])
  cv3 <- realsd[dsno]
  cv4 <- unname(perm[4][[1]])
  cv5 <- unname(perm[2][[1]])
  cv12 <- unname(perm[3][[1]])
  
  temp <- c()
  for(i in 1:20){
    temp <- c(temp,get(paste0("cv",i)))
  }
  output[k,] <- temp
}
write.csv(output,"copy.csv")
