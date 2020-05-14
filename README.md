# LATEtest
 Local Instrument Validity Tests with Causal Forests

## Description
 Assumptions sufficient for the identification of local average treatment effects (LATEs) generate necessary conditions 
    that allow to refute instrument validity. The degree of violations of instrument validity likely varies across subpopulations.
 
 `LATEtest` uses causal forests to search and test for local violations of the LATE assumptions in a data-driven way.

## Installing `LATEtest`
 To install this package in R, run the following commands:
 
  ```{js}
  library(devtools)
  install_github("farbmacher/LATEtest")
  ```

## Example
 Let `y` be the outcome, `d` an endogenous regressor, `x` an exogenous control variable and `z1 z2 z3` a set of
    potentially exogenous instruments, then the adaptive Lasso regression would be
        
    ```ruby
    n = 3000; p = 3
    cov <- matrix(c(1, 0.3, 0.3, 1), 2, 2)
    errors <-(rmvnorm(n, rep(0,2), cov))
    ```

## Author
 This package is written by Helmut Farbmacher ()

## References
 * Farbmacher H., Guber R., Klaassen S. (2020): Instrument Validity Tests with Causal Forests, working paper.
