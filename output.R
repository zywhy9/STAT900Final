library(readxl)
library(xtable)

## Plot the prior and posterior density
plotdens(n,y[1,],a,b,beta)



## Read data
data1 <- read_excel("result.xlsx", sheet = "2x1000")[,-1]
data2 <- read_excel("result.xlsx", sheet = "4x400")[,-1]
data3 <- read_excel("result.xlsx", sheet = "4x500")[,-1]
data4 <- read_excel("result.xlsx", sheet = "6x100")[,-1]
data5 <- read_excel("result.xlsx", sheet = "6x200")[,-1]
data <- rbind(data1,data2,data3,data4,data5)

## Bias vs Gelman-Rubin
plot(data$`average bias`,data$rhat1,xlab="Bias",ylab="CSRF")

## Bias vs Interval
plot(data$`average bias`,data$rhat2,xlab="Bias",ylab="Interval-based PSRF")

## MSE vs Gelman-Rubin
plot(data$`average mse`,data$rhat1,xlab="MSE",ylab="CSRF")

## MSE vs Interval
plot(data$`average mse`,data$rhat2,xlab="MSE",ylab="Interval-based PSRF")

## Data simulated
coln <- c()
for(i in 1:10){
  coln <- c(coln,paste0("Dataset ",i))
}
rownames(y) <- coln
print(xtable(y, label="table_data", digits=3, caption="Data values for analysis", align="ccccccccccc"), type="latex")

## Posterior Inference
total <- rbind(data1,data2,data3,data4,data5)
pos <- total[,c(1:5,12)]
pos <- cbind(as.vector(sapply(1:5,function(x)rep(x,10))),rep(1:10,5),pos)
colnames(pos) <- c("Design","Dataset","Real Posterior Mean", "Average Estimated Posterior Mean", 
                   "Real Posterior SD", "Average Estimated Posterior SD", 
                   "Average Bias", "Average MSE")
print(xtable(pos[1:6], label="table_pos", digits=3, caption="Posterior Inference of ", align="ccccccc"),
      type="latex", include.rownames=FALSE)

## Bias
biast <- total[,c(5:11)]
biast <- cbind(as.vector(sapply(1:5,function(x)rep(x,10))),rep(1:10,5),biast)
colnames(biast) <- c("Design","Dataset","Average Bias", "WCMCSE", 
                   "BCMCSE", "TCMCSE", "WSMCSE","BSMCSE","TSMCSE")
print(xtable(biast, label="table_bias", digits=5, caption="Summary of Posterior Bias", align="cccccccccc"),
      type="latex", include.rownames=FALSE)

## MSE
mset <- total[,c(12:18)]
mset <- cbind(as.vector(sapply(1:5,function(x)rep(x,10))),rep(1:10,5),mset)
colnames(mset) <- c("Design","Dataswet","Average MSE", "WCMCSE", 
                     "BCMCSE", "TCMCSE", "WSMCSE","BSMCSE","TSMCSE")
print(xtable(mset, label="table_mse", digits=5, caption="Summary of Posterior MSE", align="cccccccccc"),
      type="latex", include.rownames=FALSE)

## Convergence Diagnostics
diagt <- total[,19:20]
diagt <- cbind(as.vector(sapply(1:5,function(x)rep(x,10))),rep(1:10,5),diagt)
colnames(diagt) <- c("Design","Dataset","CPSF", "IPSF")
print(xtable(diagt, label="table_diag", digits=4, caption="Summary of Convergence Diagnostics", align="ccccc"),
      type="latex", include.rownames=FALSE)
