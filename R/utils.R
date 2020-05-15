#' Converts row number from tree$frame into leaf number
#'
#' @param row row
#' @param rp rp
#' @export row2leaf
row2leaf <- function(row, rp) {
  as.numeric(rownames(rp$frame)[row])
}

#' Gives the number of the parent leaf
#'
#' @param y y
#' @export parent
parent <- function(y) {
  if (y[1] != 1)
    c(Recall(if (y %% 2 == 0L) y / 2 else (y - 1) / 2), y) else y
}

#' Calculates the Group Average Treatment Effects for all nodes (terminal and non-terminal)
#'
#' @param sample data to use
#' @param treeinuse tree to use
#' @param ordered_frame ordered tree frame
#' @export gates_nodes
gates_nodes <- function(sample, treeinuse, ordered_frame) {
  leaves <- row2leaf(treeClust::rpart.predict.leaves(treeinuse, sample, type = "where"), treeinuse)
  allleaves <- lapply(leaves, parent)
  size <- nrow(treeinuse$frame)
  leafdummy <- matrix(NA, nrow = nrow(sample), ncol = size)
  nnode <- zeta <- ate_mean <- ate_se <- matrix(NA, nrow = 1, ncol = size)
  #---------------------
  for (z in 1:size) {
    l <- sort(as.numeric(rownames(treeinuse$frame)))[z]
    leafdummy[, z] <- ifelse(lapply(lapply(allleaves, is.element, l), max) == 1, 1, 0)
    subgroup <- which(as.logical(leafdummy[, z]))
    ate_mean[z] <- mean(sample$scores[subgroup])
    nnode[z] <- length(sample$scores[subgroup])
    if (sd(sample$scores[subgroup]) != 0) {
      ate_se[z] <- sd(sample$scores[subgroup]) / sqrt(nnode[z] - 1)
    } else {
      ate_se[z] <- 1
    }
    zeta[z] <- ate_mean[z] / ate_se[z]
  }
  return(list("zeta" = zeta, "ate_mean" = ate_mean, "nnode" = nnode))
}

#' Calculates the Group Average Treatment Effects for leaves only
#'
#' @param sample data to use
#' @param treeinuse tree to use
#' @param ordered_frame ordered tree frame
#' @export gates
gates <- function(sample, treeinuse, ordered_frame) {
  leaves <- row2leaf(treeClust::rpart.predict.leaves(treeinuse, sample, type = "where"), treeinuse)
  nuleaves <- length(unique(leaves))
  allleaves <- lapply(leaves, parent)
  size <- nrow(treeinuse$frame)
  leafdummy <- matrix(NA, nrow = nrow(sample), ncol = length(unique(leaves)))
  nnode <- zeta <- ate_mean <- ate_se <- matrix(NA, nrow = 1, ncol = nuleaves)
  #---------------------
  countleaf <- 0
  for (z in 1:size) {
    isleaf <- (rownames(ordered_frame)[z] %in% leaves)
    if (isleaf == TRUE) {
      countleaf <- countleaf + 1
      l <- sort(as.numeric(rownames(treeinuse$frame)))[z]
      leafdummy[, countleaf] <- ifelse(lapply(lapply(allleaves, is.element, l), max) == 1, 1, 0)
      subgroup <- which(as.logical(leafdummy[, countleaf]))
      ate_mean[countleaf] <- mean(sample$scores[subgroup])
      nnode[countleaf] <- length(sample$scores[subgroup])
      if (sd(sample$scores[subgroup]) != 0) {
        ate_se[countleaf] <- sd(sample$scores[subgroup]) / sqrt(nnode[countleaf] - 1)
      } else {
        ate_se[countleaf] <- 1
      }
      zeta[countleaf] <- ate_mean[countleaf] / ate_se[countleaf]
    }
  }
  return(list("zeta" = zeta, "ate_mean" = ate_mean, "nnode" = nnode))
}

#' Constructs pseudo-outcome, runs regression and causal forests and predicts out-of-bag (oob)
#'
#' @param est estimation sample
#' @param y relevant interval of Y
#' @param huge if set to TRUE, model for orthogonalization learned on random subset of data (size defined in tree_fraction)
#' @param tree_fraction fraction of the data used to build each tree of causal forest (and if huge==T, also used for regression forests); default=0.5
#' @param minsize causal forest insists on at least "minsize" treated and "minsize" control observations per leaf
#' @param type indicates whether we consider D=1 or D=0
#' @return pseudo-outcome, oob predictions of pseudo-outcome, oob estimates of causal forest
#' @export oob
oob <- function(est, y, huge, tree_fraction, minsize, type) {
  if (type == 1) {
    Q <- est$YinJ[, y] * est$D
  } else if (type == 0) {
    Q <- est$YinJ[, y] * (1 - est$D)
  }
  if (huge == TRUE) {
    Q.forest <- grf::regression_forest(est$X[est$huge < tree_fraction, ], Q[est$huge < tree_fraction])
    Qhat <- predict(Q.forest, est$X)$predictions
    cf <- grf::causal_forest(est$X, Q, est$Z, Y.hat = Qhat, W.hat = est$Zhat, honesty = TRUE, sample.fraction = tree_fraction, num.trees = 2000, min.node.size = minsize)
    tauhat <- predict(cf)$predictions
  } else if(huge == FALSE) {
    Qhat <- grf::regression_forest(est$X, Q)$predictions
    cf <- grf::causal_forest(est$X, Q, est$Z, Y.hat = Qhat, W.hat = est$Zhat, honesty = TRUE, num.trees = 2000, min.node.size = minsize)
    tauhat <- predict(cf)$predictions
  }
  if (type == 1) {
    Q <- -Q; Qhat <- -Qhat
  }
  return(list("Q" = Q, "Qhat" = Qhat, "tauhat" = tauhat))
}

#' Cross-fitting estimation of the GATEs
#'
#' Grow pruned tree to identify heterogeneous subgroups and picks promising leaves (based on training data)
#' Calculate magnitude of potential violations within promising leaves (in estimation sample)
#'
#' @param tra training sample
#' @param est estimation sample
#' @param y relevant interval of Y
#' @param covars provides names of covariables
#' @param type indicates whether we consider D=1 or D=0
#' @param minsize pruned tree insists on at least 2*"minsize" observations in the additional trees to identify the promising subgroups
#' @param cp sets complexity parameter which rpart uses to fit the tree before pruning; default=0
#' @return augmented inverse probability weighting scores, indicator for promising leaves
#' @export estimation
estimation <- function(tra, est, y, covars, type, minsize, cp) {
  tra$scores <- tra$tauhat + tra$Z / tra$Zhat * (tra$Q - tra$Qhat - (1 - tra$Zhat) * tra$tauhat) -
                            (1 - tra$Z) / (1 - tra$Zhat) * (tra$Q - tra$Qhat + tra$Zhat * tra$tauhat)
  formula <- paste("scores~", paste0(covars, collapse = "+"))
  tree <- rpart::rpart(formula, data = tra, method = "anova", control = rpart::rpart.control(xval = 10, minbucket = (2 * minsize), cp = cp))
  opcp <- tree$cptable[, 1][which.min(tree$cptable[, 4])]
  tree <- rpart::prune(tree, opcp)
  ordered_frame <- tree$frame[order(as.numeric(row.names(tree$frame))), ]

  estimates_tra <- gates_nodes(tra, tree, ordered_frame)
  leaves <- row2leaf(treeClust::rpart.predict.leaves(tree, tra, type = "where"), tree)
  nuleaves <- length(unique(leaves))
  uleaves <- unique(leaves)
  size <- nrow(tree$frame)
  isleaf <- relevant <- c()
  #---------------------
  if (size == 1) {
    largest_split <- TRUE
  } else {
    largest_split <- rep(FALSE, size)
    for (j in seq(2, size, 2)) {
      ifelse(estimates_tra$zeta[j] > estimates_tra$zeta[(j + 1)], largest_split[j] <- TRUE, largest_split[j + 1] <- TRUE)
    }
  }
  for (z in 1:size) {
    isleaf <- c(isleaf, (rownames(ordered_frame)[z] %in% uleaves))
    relevant <- c(relevant, (estimates_tra$zeta[z] > -qnorm(1 - 0.05 / nuleaves)))
    if (is.na(estimates_tra$zeta[z]) == TRUE) {
      relevant[z] <- FALSE
    }
  }
  estimates_tra$ate_mean <- estimates_tra$ate_mean[isleaf]
  estimates_tra$zeta <- estimates_tra$zeta[isleaf]
  estimates_tra$nnode <- estimates_tra$nnode[isleaf]

  est$scores <- est$tauhat + est$Z / est$Zhat * (est$Q - est$Qhat - (1 - est$Zhat) * est$tauhat) -
                            (1 - est$Z) / (1 - est$Zhat) * (est$Q - est$Qhat + est$Zhat * est$tauhat)
  estimates_est <- gates(est, tree, ordered_frame)
  return(list("tree" = tree, "est" = estimates_est, "tra" = estimates_tra, "relevant" = relevant[isleaf], "isleaf" = isleaf[isleaf], "largest" = largest_split[isleaf], "ordered_leaves" = ordered_frame[isleaf, ]))
}
