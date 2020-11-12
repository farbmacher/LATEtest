# LATEtest
Local Instrument Validity Tests with Causal Forests

## Description
Assumptions sufficient for the identification of local average treatment effects (LATEs) generate necessary conditions 
that allow to refute instrument validity. The degree of violations of instrument validity likely varies across 
subpopulations.
 
`LATEtest` uses causal forests to search and test for local violations of the LATE assumptions in a data-driven way.

## Installing `LATEtest`
To install this package in R, run the following commands:

```R
install.packages("devtools")
library(devtools)
install_github("farbmacher/LATEtest")
```

## Example
Let `Y` be the outcome, `D` an endogenous treatment indicator, `Z` a binary instrument, `X` predetermined covariates:
        
```R
library("LATEtest")
library("rpart.plot")

n = 3000; p = 3; rho=0.3
u <- rnorm(n)
v <- rho*u + sqrt(1-rho^2)*rnorm(n)
X <- matrix(rnorm(n * p), n, p)
colnames(X) <- paste("Xvar", 1:p, sep="")
Z <- rbinom(n, size = 1, prob = 0.5)
D<-as.numeric(0.2 * Z + v > 0)

# Local violation of the exclusion restriction:
gamma <- as.numeric(ifelse(X[, 2] < -1, 1.25, 0))
Y <- as.vector(D + gamma * Z + u)
data <- as.data.frame(cbind(Y,D,Z,X))

# Perform test:
covars = paste0(colnames(data)[4:ncol(data)])
test <- LATEtest(data = data, covars = covars, subsets = 4, alpha = 0.05)
test

# Draw plot of pruned tree that led to local violation of LATE assumptions:
maxtree <- eval(parse(text = paste("test$treelist$tree_", test$maxTtree$label, test$maxTtree$J,sep = "")))
rpart.plot(maxtree, extra = 101, box.palette = "GyRd", shadow.col = "gray", nn = TRUE, roundint = FALSE)
```

## Author
This package is written by Helmut Farbmacher (farbmacher@econ.lmu.de). I thank Christina Niessl for excellent research assistance.

## References
* Farbmacher H., Guber R., Klaassen S. (2020): Instrument Validity Tests with Causal Forests, *Journal of Business and Economic Statistics*, [Link](https://www.tandfonline.com/doi/abs/10.1080/07350015.2020.1847122?journalCode=ubes20).
