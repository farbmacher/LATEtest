#' LATEtest
#'
#' Main function: Prepares data, calls help functions, stores results.
#' See the LATEtest website on \url{https://github.com/farbmacher/LATEtest} for more information, documentation and examples.
#' Please report any bugs to \email{farbmacher@@econ.lmu.de}
#'
#' @param data provides dataset: outcome must be labelled Y, treatment D and instrument Z
#' @param covars provides names of covariables in data (e.g., "Xvar1" "Xvar2" "Xvar3")
#' @param huge if set to TRUE, model for orthogonalization learned on random subset of data (size defined in tree_fraction)
#' @param tree_fraction fraction of the data used to build each tree of causal forest (and if huge==T, also used for regression forests); default=0.5
#' @param minsize causal forest insists on at least "minsize" treated and "minsize" control observations per leaf, and pruned tree insists on at least 2*"minsize" observations in the additional trees to identify the promising subgroups
#' @param cp sets complexity parameter which rpart uses to fit the tree before pruning; default=0
#' @param subsets number of subsets used to discretize the outcome (using an equidistant grid from the minimum to the maximum value of Y)
#' @param alpha nominal significance level of the test
#' @param seed set random seed number
#' @return list of the pruned trees, test results
#' @examples
#' \dontrun{
#' library("LATEtest")
#' library("mvtnorm", "rpart.plot")
#'
#' # Generate data:
#' n = 3000; p = 3
#' cov <- matrix(c(1, 0.3, 0.3, 1), 2, 2)
#' errors <-(rmvnorm(n, rep(0,2), cov))
#' X <- matrix(rnorm(n * p), n, p)
#' colnames(X) <- paste("Xvar", 1:p, sep="")
#' Z <- rbinom(n, size = 1, prob = 0.5)
#' D<-as.numeric(0.2 * Z + errors[, 1] > 0)
#' # local violation of the exclusion restriction:
#' gamma <- as.numeric(ifelse(X[, 2] < -1, 1.25, 0))
#' Y <- as.vector(D + gamma * Z + errors[, 2])
#' data <- as.data.frame(cbind(Y,D,Z,X))
#'
#' # Perform test:
#' covars = paste0(colnames(data)[4:ncol(data)])
#' test <- LATEtest(data = data, covars = covars, subsets = 4, alpha = 0.05)
#' test
#'
#' # Draw plot of pruned tree that led to local violation of LATE assumptions:
#' maxtree <- eval(parse(text = paste("test$treelist$tree_", test$maxTtree$label, test$maxTtree$J,sep = "")))
#' rpart.plot::rpart.plot(maxtree, extra = 101, box.palette = "GyRd", shadow.col = "gray", nn = TRUE, roundint = FALSE)
#' }

LATEtest <- function(data, covars, huge=FALSE, tree_fraction=0.5, minsize=100, cp=0, subsets=4, alpha=0.05, seed=10101) {
  set.seed(seed)

  # Prepare data
  ##############
  n <- nrow(data)
  Yorig <- data$Y
  Y <- as.numeric(cut(Yorig, breaks = seq(from = min(Yorig) - 0.001, to = max(Yorig) + 0.001, length.out = subsets + 1)))
  luY <- length(unique(Y))
  YinJ <- matrix(NA, nrow = n, ncol = luY)
  for (y in 1:luY) {
    YinJ[, y] <- ifelse(Y == sort(unique(Y))[y], 1, 0)
  }
  data$YinJ <- YinJ
  X <- as.matrix(data[, covars])
  data$X <- X
  if (huge == TRUE) {
    data$huge <- runif(n)
  }
  trID <- which(data$Z == 1); conID <- which(data$Z == 0)
  index <- c(sample(trID, length(trID) * 0.5), sample(conID, length(conID) * 0.5))
  sampleA <- data[index, ]; sampleB <- data[-index, ]
  if (huge == TRUE) {
    Z.forest <- grf::regression_forest(sampleA$X[sampleA$huge < tree_fraction, ], sampleA$Z[sampleA$huge < tree_fraction])
    sampleB$Zhat <- predict(Z.forest, sampleB$X)$predictions
    Z.forest <- grf::regression_forest(sampleB$X[sampleB$huge < tree_fraction, ], sampleB$Z[sampleB$huge < tree_fraction])
    sampleA$Zhat <- predict(Z.forest, sampleA$X)$predictions
  } else if (huge == FALSE) {
    Z.forest <- grf::regression_forest(sampleA$X, sampleA$Z); sampleB$Zhat <- predict(Z.forest, sampleB$X)$predictions
    Z.forest <- grf::regression_forest(sampleB$X, sampleB$Z); sampleA$Zhat <- predict(Z.forest, sampleA$X)$predictions
  }

  # Prepare results
  #################
  zeta_all <- sel1_all <- sel2_all <- c()
  leafinfo <- treelist <- vector("list", length = 1)
  zeta <- sel1 <- sel2 <- treelist_names <- c()
  zetaQ1_info <- vector("list", length = (2 * luY)); zetaQ0_info <- vector("list", length = (2 * luY))
  treelistQ1 <- vector("list", length = luY); treelistQ0 <- vector("list", length = luY)

  #------------------------------------------------------------------------------------------------
  # first for Q1 (d=1):
  for (y in 1:luY) {

    temp_oob <- oob(sampleA, y, huge, tree_fraction, minsize, 1)
    sampleA$Q <- temp_oob$Q; sampleA$Qhat <- temp_oob$Qhat; sampleA$tauhat <- temp_oob$tauhat
    rm(temp_oob)
    temp_oob <- oob(sampleB, y, huge, tree_fraction, minsize, 1)
    sampleB$Q <- temp_oob$Q; sampleB$Qhat <- temp_oob$Qhat; sampleB$tauhat <- temp_oob$tauhat
    rm(temp_oob)

    # Sample A using treeB #
    ########################
    temp <- estimation(sampleB, sampleA, y, covars, 1, minsize, cp)
    zeta <- cbind(zeta, temp$est$zeta)
    sel1 <- c(sel1, temp$relevant)                    # exclude leaves with a significant fraction of compliers
    sel2 <- c(sel2, temp$relevant * temp$largest)     # take only promising leaves (largest binary split)

    zetaQ1_info[[(2*y-1)]] <- data.frame(label="Q1A", J=y, knot=rownames(temp$ordered_leaves), n_train=temp$tra$nnode,
                                         Tstat_train=temp$tra$zeta, is_relevant=temp$relevant,
                                         is_largest=temp$largest, sep="| ", n_est=t(temp$est$nnode), Tstat_est=t(temp$est$zeta))
    treelistQ1[[(2*y-1)]] <- temp$tree
    treelist_names <- c(treelist_names, paste0("tree_Q1A", y))
    rm(temp)

    # Sample B using treeA #
    ########################
    temp <- estimation(sampleA, sampleB, y, covars, 1, minsize, cp)
    zeta <- cbind(zeta, temp$est$zeta)
    sel1 <- c(sel1, temp$relevant)                    # exclude leaves with a significant fraction of compliers
    sel2 <- c(sel2, temp$relevant * temp$largest)     # take only promising leaves (largest binary split)

    zetaQ1_info[[(2*y)]] <- data.frame(label="Q1B", J=y, knot=rownames(temp$ordered_leaves), n_train=temp$tra$nnode,
                                       Tstat_train=temp$tra$zeta, is_relevant=temp$relevant,
                                       is_largest=temp$largest, sep="| ", n_est=t(temp$est$nnode), Tstat_est=t(temp$est$zeta))
    treelistQ1[[(2*y)]] <- temp$tree
    treelist_names <- c(treelist_names, paste0("tree_Q1B", y))
    rm(temp)

  }

  #------------------------------------------------------------------------------------------------
  # now for Q0 (d=0):
  for (y in 1:luY) {

    temp_oob <- oob(sampleA, y, huge, tree_fraction, minsize, 0)
    sampleA$Q <- temp_oob$Q; sampleA$Qhat <- temp_oob$Qhat; sampleA$tauhat <- temp_oob$tauhat
    rm(temp_oob)
    temp_oob <- oob(sampleB, y, huge, tree_fraction, minsize, 0)
    sampleB$Q <- temp_oob$Q; sampleB$Qhat <- temp_oob$Qhat; sampleB$tauhat <- temp_oob$tauhat
    rm(temp_oob)

    # Sample A using treeB #
    ########################
    temp <- estimation(sampleB, sampleA, y, covars, 0, minsize, cp)
    zeta <- cbind(zeta, temp$est$zeta)
    sel1 <- c(sel1, temp$relevant)                    # exclude leaves with a significant fraction of compliers
    sel2 <- c(sel2, temp$relevant * temp$largest)     # take only promising leaves (largest binary split)

    zetaQ0_info[[(2*y-1)]] <- data.frame(label="Q0A", J=y, knot=rownames(temp$ordered_leaves), n_train=temp$tra$nnode,
                                         Tstat_train=temp$tra$zeta, is_relevant=temp$relevant,
                                         is_largest=temp$largest, sep="| ", n_est=t(temp$est$nnode), Tstat_est=t(temp$est$zeta))
    treelistQ0[[(2*y-1)]] <- temp$tree
    treelist_names <- c(treelist_names, paste0("tree_Q0A", y))
    rm(temp)

    # Sample B using treeA #
    ########################

    temp <- estimation(sampleA, sampleB, y, covars, 0, minsize, cp)
    zeta <- cbind(zeta, temp$est$zeta)
    sel1 <- c(sel1, temp$relevant)                    # exclude leaves with a significant fraction of compliers
    sel2 <- c(sel2, temp$relevant * temp$largest)     # take only promising leaves (largest binary split)

    zetaQ0_info[[(2*y)]] <- data.frame(label="Q0B", J=y, knot=rownames(temp$ordered_leaves), n_train=temp$tra$nnode,
                                       Tstat_train=temp$tra$zeta, is_relevant=temp$relevant,
                                       is_largest=temp$largest, sep="| ", n_est=t(temp$est$nnode), Tstat_est=t(temp$est$zeta))
    treelistQ0[[(2*y)]] <- temp$tree
    treelist_names <- c(treelist_names, paste0("tree_Q0B", y))
    rm(temp)

  }
  zeta_all <- cbind(zeta_all, zeta); sel1_all <- as.logical(t(cbind(sel1_all, sel1))); sel2_all <- as.logical(t(cbind(sel2_all, sel2)))

  #-------------------------------------------------------------------------------------------
  # Store tree info :
  treelist <- c(treelistQ1, treelistQ0)
  names(treelist) <- treelist_names
  zetaQ1_info <- do.call("rbind", zetaQ1_info); zetaQ0_info <- do.call("rbind", zetaQ0_info)
  leafinfo <- rbind(zetaQ1_info, zetaQ0_info)
  rownames(leafinfo) <- paste(1:nrow(leafinfo))
  #-------------------------------------------------------------------------------------------

  #------------------------------------------------------------------------------------------------
  ## test procedure
  #------------------------------------------------------------------------------------------------
  sel_na <-  c()
  for (x in 1:ncol(zeta_all)) {
    if (is.na(zeta_all[x]) == F) {
      sel_na <- c(sel_na, x)
    }
  }
  zeta <- t(zeta_all[sel_na]); sel1_all <- t(sel1_all[sel_na]); sel2_all <- t(sel2_all[sel_na])

  p0 <- ncol(zeta)
  T0 <- max(zeta)
  cv0 <- qnorm(1 - alpha / p0)

  sel1_final <- c()
  for (x in 1:ncol(zeta)) {
    if (sel1_all[x] == 1) {
      sel1_final <- c(sel1_final, x)
    }
  }
  zeta1 <- t(zeta[sel1_final])
  p1 <- ncol(zeta1)
  T1 <- max(zeta1)
  cv1 <- qnorm(1 - alpha / p1)

  sel2_final <- c()
  for(x in 1:ncol(zeta)) {
    if (sel2_all[x] == 1) {
      sel2_final <- c(sel2_final, x)
    }
  }
  zeta2 <- t(zeta[sel2_final])
  p2 <- ncol(zeta2)
  T2 <- max(zeta2)
  cv2 <- qnorm(1 - alpha / p2)

  #-------------------------------------------------------------------------------------------
  # Store additional tree info:
  maxTtree <- c()
  leafinfo$Tmax <- ifelse(leafinfo$Tstat_est == T2, 1, 0)
  if (mean(leafinfo$Tmax, na.rm = T) != 0) {
    ident <- leafinfo$Tmax == 1
    ident <- ifelse(is.na(ident), FALSE, ident)
    maxTtree$label <- leafinfo$label[ident]
    maxTtree$J <- leafinfo$J[ident]
  }
  #-------------------------------------------------------------------------------------------
  # Store chosen parameters:
  parameters <- c()
  parameters$covars <- covars; parameters$subsets <- subsets; parameters$minsize <- minsize
  parameters$alpha <- alpha; parameters$huge <- huge
  if (huge == TRUE) {
    parameters$tree_fraction <- tree_fraction
  }
  #-------------------------------------------------------------------------------------------
  descriptives <- c()
  descriptives$meanY_D_Z <- c(summary(Yorig)[4], summary(data$D)[4], summary(data$Z)[4])
  descriptives$min_mean_maxZhatA <- c(summary(sampleA$Zhat)[1], summary(sampleA$Zhat)[4], summary(sampleA$Zhat)[6])
  descriptives$min_mean_maxZhatB <- c(summary(sampleB$Zhat)[1], summary(sampleB$Zhat)[4], summary(sampleB$Zhat)[6])
  results <- c()
  results$nu_ineq <- cbind(p0, p1, p2)
  results$Tmax <- cbind(T0, T1, T2)
  results$cv <- cbind(cv0, cv1, cv2)
  results$reject <- cbind(T0 > cv0, T1 > cv1, T2 > cv2)
  results$pvalue <- cbind(min(p0 * (1 - pnorm(T0)), 1), min(p1 * (1 - pnorm(T1)), 1), min(p2 * (1 - pnorm(T2)), 1))
  return(list("treelist" = treelist, "maxTtree" = maxTtree, "descriptives" = descriptives, "parameters" = parameters,
              "leafinfo" = leafinfo, "results" = results))
}
