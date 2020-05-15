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
  data <- as.data.frame(cbind(Y,D,Z,X))   # if we add further variables to data, then the call of LATEtest.R in the simulation file has to be adapted!
  return(data)
}

