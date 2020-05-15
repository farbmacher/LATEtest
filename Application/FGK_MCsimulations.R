#################################################
# Simulations from
# "Instrument Validity Tests with Causal Forests"
# Helmut Farbmacher et al 2020
# May 2020
#################################################

rm(list=ls())
packages <- c("LATEtest","foreach", "doParallel", "mvtnorm")
invisible(lapply(packages, library, character.only = TRUE))

setwd("/Users/helmut/Desktop/Simulations")   # set working directory path

######################################################################################
# Function for data simulation
fct_datasim <- function(setup,dgp){
  p <- 3
  betaXY <- c(0.3,0.3,0.3)
  cov <- matrix(c(1,0.3,0.3,1),2,2)
  errors <-(rmvnorm(n,rep(0,2),cov))
  X <- matrix(rnorm(n*p), n, p)
  colnames(X) <- paste("Xvar", 1:p, sep="")

  if(dgp==0) {          # no violations of A1 & A3; inequalities are binding
    a <- 0; alpha =  rep(a, n)
    g <- 0; gamma = rep(g, n)
  } else if(dgp==1) {   # no violations of A1 & A3 but with heterogeneity; inequaltities are not binding (all tests are conservative)
    a <- 0.2; alpha =  rep(a, n)
    g <- 0;   gamma = rep(g, n)
  } else if(dgp==2) {   # local violation of monotonicity
    alpha <- as.numeric(ifelse(X[,1]< qnorm(0.15), -0.75, 0.4))
    g <- 0; gamma = rep(g, n)
  } else if(dgp==3) {   # local violation of exclusion restriction
    a <- 0.2; alpha =  rep(a, n)
    gamma <- as.numeric(ifelse(X[,2]< qnorm(0.15), 1.25, 0))
  } else if(dgp==4) {   # global violation of exclusion restriction
    a <- 0.2; alpha =  rep(a, n)
    g <- 0.5; gamma = rep(g, n)
  } else if(dgp==5) {   # global violation of exclusion restriction but with sign heterogeneity
    a <- 0.2; alpha =  rep(a, n)
    gamma <- as.numeric(ifelse(X[,3]< qnorm(0.5), 0.5, -0.5))
  } else {
    stop("invalid choice of dgp")
  }

  if (setup == 'A'){         # basic setup (randomized experiment)
    Z <- rbinom(n, size=1, prob=0.5)
    Dlatent <- alpha*Z+errors[,1]
    D<-as.numeric( Dlatent>0 )
  } else if (setup == 'B'){  # easy propensity and strong confounding
    Zlatent = 0.2*X[,1] + 0.2*X[,2] + 0.2*X[,3] + rnorm(n)
    Z = as.numeric( Zlatent>0 )
    b = 0.5 * log(1 + exp(X[,1] + X[,2] + X[,3]))
    Dlatent = b + alpha*Z + errors[,1]
    D = as.numeric( Dlatent > 0 )
  } else {
    stop("invalid choice of setup")
  }
  Y <- as.vector(1*D+gamma*Z+X%*%betaXY+errors[,2])
  data <- as.data.frame(cbind(Y,D,Z,X))
  return(data)
}
######################################################################################

seed=1234
set.seed(seed)

##### Setup DGP #####
#####################
setup <- "A"          # A: randomized experiment, B: easy confounding of Z and strong confounding of D
n <- 3000             # number of observations
R <- 100              # number of Monte Carlo replications

subsets <- 4
siglevel <- 0.05

##### Setup parallel computing #####
####################################
system <- Sys.info()['sysname']
if(system == 'Windows') {
  cl <- makeCluster(detectCores())
  clusterEvalQ(cl, c(library(LATEtest),library(grf),library(rpart),library(treeClust),library(mvtnorm)))
} else if(system == 'Darwin') {
  cl <- makeForkCluster(detectCores())
}
registerDoParallel(cl)
showConnections()

sink(file=paste("FGK","_setup",setup,".txt",sep=""), split=TRUE)
print(Sys.time())
print("-----------------------")
print("setup, n, R, seed, subsets, siglevel")
print(paste(setup, n, R, seed, subsets, siglevel, sep=" / "))
print("--------------------------------------------")
print("--------------------------------------------")

##### Run simulation #####
##########################
for (dgp in 0:5) {         #loop over all dgps; see fct_datasim() above
  start_time = Sys.time()
  sim <- foreach(r=1:R, .combine=rbind) %dopar% {
    set.seed(r)
    data <- fct_datasim(setup=setup, dgp=dgp)

    test <- LATEtest(data=data, covars=paste0(colnames(data)[4:ncol(data)]), subsets=subsets, alpha=siglevel)
    return(c(r, test$results$Tmax,test$results$cv,test$results$reject,test$results$nu_ineq))
  }
  end_time = Sys.time()
  print(end_time - start_time)
  colnames(sim) <- c("round","Tmax0","Tmax1","Tmax2","cv0","cv1","cv2","reject0","reject1","reject2","dim0","dim1","dim2")
  sim <- as.data.frame(sim)

  ##### Display results #####
  ###########################
  print("dgp")
  print(dgp)
  options(digits=3); print(head(sim)); options(digits=4)
  print("summary(sim$Tmax0/1/2)")
  print(rbind(summary(sim$Tmax0),summary(sim$Tmax1),summary(sim$Tmax2)))
  print("summary(sim$cv0/1/2)")
  print(rbind(summary(sim$cv0),summary(sim$cv1),summary(sim$cv2)))
  print("mean(sim$reject0/1/2)")
  print(c(mean(sim$reject0),mean(sim$reject1),mean(sim$reject2)))
  print("")
  print("summary(sim$dim0/1/2)")
  print(rbind(summary(sim$dim0),summary(sim$dim1),summary(sim$dim2)))
  print("table(sim$dim0)")
  print(table(sim$dim0))
  print("table(sim$dim1)")
  print(table(sim$dim1))
  print("table(sim$dim2)")
  print(table(sim$dim2))
  save.image(file=paste("FGK","_setup",setup,"_dgp",dgp,".Rdata",sep=""))
  print("--------------------------------------")
}
sink(file=NULL)
stopCluster(cl)
registerDoParallel()

