# reinforced risk prediction with budget constraint: variable selection version -------------------


#' Model fit for the training set, variable selection version
#'
#' \code{modelFit_VS} outputs the FPCA (functional principal component analysis) decomposition and the elastic net logistic regression
#' at each time grid. This function is used when the baseline covariates Z are high-dimensional.
#'
#' @export
#' @param Y The outcome variable, vector of length \eqn{n}, taking values in \eqn{{1, 0, NA}}, where 1 = disease, 0 = not, NA = missing.
#' @param X Observed longitudinal biomarker, matrix of \eqn{n} by \eqn{nTotal}, where \eqn{nTotal} denotes the total number of time grids.
#' Missing values are denoted by NA.
#' @param Z Other baseline covariates.
#' @param startT Time of the first prediction, denoted by \eqn{t_1} in the manuscript. For instance, if the time grids are \eqn{{0,1/60,2/60,...,1}},
#' then startT = 25 means that the first prediction is made at \eqn{t = 24/60}.
#' @param pve Proportion of variance explained in FPCA.
#' @param nbasis Number of B-spline basis functions needed for estimation of the mean function and smoothing of covariance.
#' @param weight Weight for each individual.
#'
#' @return
#' \item{list_fpcaFit}{FPCA decomposition at each time grid from startT to the end.}
#' \item{list_cvfit}{Elastic net logistic regression at each time grid from startT to the end.}
#'
#' @examples
#' \donttest{
#' library(reinforcedPred)
#'
#' # take the example training data (high dimensional Z) from the reinforcedPred package
#' # see documentation for details about the data set train_data_mulZ
#' Y <- as.numeric(train_data_mulZ$Y)
#' tildeX.missing <- as.matrix(train_data_mulZ[,2:62])
#' Z <- as.matrix(train_data_mulZ[,63:dim(train_data_mulZ)[2]])
#'
#' # analysis starts
#' startT <- 25
#' weight <- rep(1, length(Y))
#'
#' result <- modelFit_VS(Y, tildeX.missing, Z, startT, pve = 0.99, nbasis = 10, weight)
#'
#' # obtained elastic net logistic regression fit and FPCA decompositions
#' list_cvfit <- result$list_cvfit
#' list_fpcaFit <- result$list_fpcaFit
#' }


modelFit_VS <- function(Y, X, Z, startT, pve, nbasis, weight) {

  n <- length(Y)                                 # sample size
  nTotal <- dim(X)[2]                            # number of total time grid

  list_fpcaFit <- list()                         # FPCA Fit at each time grid
  list_cvfit <- list()                           # parameter estimates at each time grid

  # curT denotes the current time grid #
  for (curT in startT:nTotal) {

    curX <- X[,1:curT]                           # the up-to-date longitudinal data
    fpcaFit <- refund::fpca.sc(Y = curX, nbasis = nbasis, var = TRUE, pve = pve)
    list_fpcaFit <- append(list_fpcaFit, list(fpcaFit))

    # number of FPCs (functional principal components) #
    K <- fpcaFit$npc

    # FPC scores #
    FPC.score <- fpcaFit$scores

    # fit an elastic net type logistic regression using FPC scores and baseline covariates as predictors #
    x <- cbind(FPC.score, Z)

    # subjects who have no missing values for the outcome and the covariates #
    cc <- stats::complete.cases(Y) & stats::complete.cases(x)

    Y.cc <- Y[cc]
    x.cc <- x[cc, ]
    weight.cc <- weight[cc]

    cvfit <- glmnet::cv.glmnet(x.cc, Y.cc, weights = weight.cc, family = "binomial", type.measure = "class")

    list_cvfit <- append(list_cvfit, list(cvfit))

  }

  return (list("list_fpcaFit" = list_fpcaFit, "list_cvfit" = list_cvfit))
}


#' Risk prediction on the test set, variable selection version
#'
#' For a fixed threshold value \eqn{\tau}, \code{modelPredict_VS} predicts the outcome \eqn{Y}  for subjects in the test set.
#' This function also outputs the cost associated with the prediction procedure. This function is used when the baseline
#' covariates Z are high-dimensional.
#'
#' @export
#' @param list_fpcaFit Obtained FPCA decomposition from \code{modelFit_VS}.
#' @param list_cvfit Obtained elastic net logistic regression from \code{modelFit_VS}.
#' @param Xtest Longitudinal biomarker data for subjects in the test set, matrix of \eqn{testn} by \eqn{nTotal}. Missing values
#' are denoted by NA.
#' @param Ztest Other baseline covariates for subjects in the test set.
#' @param startT Time of the first prediction, denoted by \eqn{t_1} in the manuscript. For instance, if the time grids are \eqn{{0,1/60,2/60,...,1}},
#' then startT = 25 means that the first prediction is made at \eqn{t = 24/60}.
#' @param tau The threshold value \eqn{\tau}.
#'
#' @return
#' \item{final.label}{Predicted outcome \eqn{Y} for subjects in the test set, vector of length \eqn{testn}.}
#' \item{avg.cost}{Average cost when we applied this prediction procedure to the test set. }
#' \item{cost}{Cost for each subject, vector of length \eqn{testn}. For some subjects, we make a definite decision early. For others, we follow up with a long period of time. Hence the cost is different for each individual. }
#'
#' @examples
#' # see the example from function reinforced_VS.


modelPredict_VS <- function(list_fpcaFit, list_cvfit, Xtest, Ztest, startT, tau) {

  testn <- dim(Xtest)[1]                         # number of subjects in the test data set
  nTotal <- dim(Xtest)[2]                        # number of total time grid

  prediction <- matrix(NA, nrow = testn, ncol = nTotal)

  # curT denotes the current time grid #
  # from starting time grid to time grid (nTotal - 1), we use classification with reject option #
  for (curT in startT:(nTotal - 1)) {

    # the FPCA fit at the current time grid #
    fpcaFit <- list_fpcaFit[[curT - startT + 1]]

    # number of FPCs #
    K <- fpcaFit$npc

    # FPC basis functions #
    FPC.basis <- fpcaFit$efunctions

    # estimated mean function #
    mu <- fpcaFit$mu

    # estimated eigenvalues of the covariance operator, i.e., variances of FPC scores #
    evalues <- fpcaFit$evalues

    # estimated measurement error variance #
    sigma2 <- fpcaFit$sigma2

    curXtest <- Xtest[,1:curT]                   # the up-to-date longitudinal measurement

    # estimated FPC scores for subjects in the test data set #
    estimated.FPC.score <- matrix(NA, nrow = testn, ncol = K)
    for (i in 1:testn) {
      ind <- which(!is.na(curXtest[i, ]))
      temp.matrix <- matrix(FPC.basis[ind, ], ncol = K)

      if (K == 1) {
        estimated.FPC.score[i, ] <- MASS::ginv(t(temp.matrix) %*% temp.matrix + sigma2 * evalues) %*% t(temp.matrix) %*% (Xtest[i, ind] - mu[ind])
      } else {
        estimated.FPC.score[i, ] <- MASS::ginv(t(temp.matrix) %*% temp.matrix + sigma2 * diag(evalues)) %*% t(temp.matrix) %*% (Xtest[i, ind] - mu[ind])
      }

    }

    # the parameter estimates at the current time grid #
    cvfit <- list_cvfit[[curT - startT + 1]]

    # covariates including the estimated FPC scores and other baseline covariates
    testx <- cbind(estimated.FPC.score, Ztest)

    # predict the label Y #
    f <- stats::predict(cvfit, newx = testx, s = "lambda.1se", type = "link")

    # 1 = disease, 0 = not, -100 = uncertain (reject option) #
    predict.Y <- rep(-100, testn)
    ind1 <- which(f > tau)
    ind0 <- which(f < -tau)
    predict.Y[ind1] <- 1
    predict.Y[ind0] <- 0

    prediction[, curT] <- predict.Y
  }

  # for the final time grid, we only use the binary classification #
  fpcaFit <- list_fpcaFit[[nTotal - startT + 1]]

  # number of FPCs #
  K <- fpcaFit$npc

  # FPC basis functions #
  FPC.basis <- fpcaFit$efunctions

  # estimated mean function #
  mu <- fpcaFit$mu

  # estimated eigenvalues of the covariance operator, i.e., variances of FPC scores #
  evalues <- fpcaFit$evalues

  # estimated measurement error variance #
  sigma2 <- fpcaFit$sigma2

  # estimated FPC scores for subjects in the test data #
  estimated.FPC.score <- matrix(NA, nrow = testn, ncol = K)
  for (i in 1:testn) {
    ind <- which(!is.na(Xtest[i, ]))
    temp.matrix <- matrix(FPC.basis[ind, ], ncol = K)

    if (K == 1) {
      estimated.FPC.score[i, ] <- MASS::ginv(t(temp.matrix) %*% temp.matrix + sigma2 * evalues) %*% t(temp.matrix) %*% (Xtest[i, ind] - mu[ind])
    } else {
      estimated.FPC.score[i, ] <- MASS::ginv(t(temp.matrix) %*% temp.matrix + sigma2 * diag(evalues)) %*% t(temp.matrix) %*% (Xtest[i, ind] - mu[ind])
    }
  }

  # the parameter estimates at the last time grid #
  cvfit <- list_cvfit[[nTotal - startT + 1]]

  # covariates including the estimated FPC scores and other baseline covariates
  testx <- cbind(estimated.FPC.score, Ztest)

  # predict the label Y #
  f <- stats::predict(cvfit, newx = testx, s = "lambda.1se", type = "link")

  prediction[, nTotal] <- ifelse(stats::pnorm(f) > 0.5, 1, 0)

  cost <- rep(NA, testn)
  final.label <- rep(NA, testn)

  for(i in 1:testn) {
    cost[i] <- suppressWarnings(min(min(which(prediction[i, ] == 1)), min(which(prediction[i, ] == 0))))
    final.label[i] <- prediction[i, cost[i]]
  }

  avg.cost <- mean(cost)

  return (list("final.label" = final.label, "avg.cost" = avg.cost, "cost" = cost))
}


#' Reinforced risk prediction with budget constraint, variable selection version
#'
#' \code{reinforced_VS} implements a cross-validation approach to find an optimal \eqn{\tau} such that the
#' misclassification error is minimized under a certain budget constraint. This function is used when the baseline
#' covariates are of high-dimension.
#' @export
#' @param Y The outcome variable, vector of length \eqn{n}, taking values in \eqn{{1, 0, NA}}, where 1 = disease, 0 = not, NA = missing.
#' @param X Observed longitudinal biomarker, matrix of \eqn{n} by \eqn{nTotal}, where \eqn{nTotal} denotes the total number of time grids.
#' Missing values are denoted by NA.
#' @param Z Other baseline covariates.
#' @param budget The budget constraint. For instance, if the time grids are \eqn{{0,1/60,2/60,...,1}}. Budget = 30 means that the average follow up was no longer than 30 time grids.
#' This is equivalent to saying that on average, we want to make a definite prediction before time \eqn{t = 0.5}.
#' @param folds Folds in cross-validation, usually 5 or 10.
#' @param startT Time of the first prediction, denoted by \eqn{t_1} in the manuscript. For instance, if the time grids are \eqn{{0,1/60,2/60,...,1}},
#' then startT = 25 means that the first prediction is made at \eqn{t = 24/60}.
#' @param pve Proportion of variance explained in FPCA, default value is 0.99.
#' @param nbasis Number of B-spline basis functions needed for estimation of the mean function and smoothing of covariance.
#' Default value is 10 in refund package, sometimes a smaller number is needed when there are a small number of time grids.
#' @param weight A user-supplied weight for each individual. If the user did not supply the weight, we use an inverse probability
#' weighting method to calculate a weight. See details in section 3.4 of the manuscript.
#'
#' @return
#' \item{final.result}{The FPCA fit and the elastic net logistic regression fit at each time grid from startT to the end.}
#' \item{final.tau}{The optimal \eqn{\tau} that minimizes the misclassification error under the budget constraint.}
#'
#' @examples
#' \donttest{
#' library(reinforcedPred)
#' set.seed(1)
#'
#' # take the example training data (high dimensional Z) from the reinforcedPred package
#' # see documentation for details about the data set train_data_mulZ
#' Y <- as.numeric(train_data_mulZ$Y)
#' tildeX.missing <- as.matrix(train_data_mulZ[,2:62])
#' Z <- as.matrix(train_data_mulZ[,63:dim(train_data_mulZ)[2]])
#'
#' # analysis starts
#' budget <- 45
#' folds <- 5
#' startT <- 25
#'
#' result <- reinforced_VS(Y, tildeX.missing, Z, budget, folds, startT, pve = 0.99, nbasis = 10)
#'
#' # obtained elastic net logistic regression fit and FPCA decompositions
#' list_cvfit <- (result$final.result)$list_cvfit
#' list_fpcaFit <- (result$final.result)$list_fpcaFit
#'
#' # optimal tau that minimizes the misclassification error under the budget constraint
#' final.tau <- result$final.tau
#' final.tau
#'
#' # use the fitted model to predict the label Y for subjects in the test data
#' # see documentation for details about the data set test_data_mulZ
#' testY <- as.numeric(test_data_mulZ$testY)
#' test.tildeX.missing <- as.matrix(test_data_mulZ[,2:62])
#' test.Z <- as.matrix(test_data_mulZ[,63:dim(test_data_mulZ)[2]])
#'
#' pred <- modelPredict_VS(list_fpcaFit, list_cvfit, test.tildeX.missing, test.Z, startT, final.tau)
#'
#' # predicted outcome Y for each subject in the test data
#' predY.test <- pred$final.label
#' # misclassification error
#' mis.error <- sum(predY.test != testY, na.rm = TRUE) / sum(!is.na(testY))
#' mis.error
#'
#' # the average cost when we applied the prediction procedure to the test data
#' pred$avg.cost
#' }


reinforced_VS <- function(Y, X, Z, budget, folds, startT, pve = 0.99, nbasis = 10, weight) {

  n <- length(Y)                  # sample size

  # if the user does not supply the weight #
  if (nargs() < 9) {

    # FPCA based on all longitudinal information #
    full.fpcaFit <- refund::fpca.sc(Y = X, var = TRUE, pve = pve)

    # FPC scores: this term will be later used to model the probability that Y is missing #
    full.FPC.score <- full.fpcaFit$scores

    R <- !is.na(Y)                 # indicator variable for the non-missingness of the outcome variable Y #

    if (sum(R) == length(Y)) {
      # if there is no missing value for Y #
      weight <- rep(1, n)
    } else {
      # if outcome is missing for some subjects #

      # model the probability that the outcome is missing #
      fit.propen <- stats::glm(R ~ full.FPC.score + Z, family = "binomial")
      prob <- fit.propen$fitted.values
      weight <- 1 / prob
    }

  }

  # we search over a sequence of tau values #
  tau.sequence <- seq(0, 3, 0.02)

  # record the misclassification error and the cost, each row is for fixed fold with different taus #
  error.matrix <- matrix(NA, nrow = folds, ncol = length(tau.sequence))
  cost.matrix <- matrix(NA, nrow = folds, ncol = length(tau.sequence))

  # cross-validation procedure #
  samplen = sample(n)

  for (i in 1:folds) {
    tst_idx = samplen[seq(i, n, by = folds)]                   # index set for test data
    trn_idx = setdiff(1:n, tst_idx)                            # index set for training data

    trn_Y <- Y[trn_idx]
    trn_X <- X[trn_idx, ]
    trn_weight <- weight[trn_idx]

    tst_Y <- Y[tst_idx]
    tst_X <- X[tst_idx, ]

    if (is.vector(Z)) {
      trn_Z <- Z[trn_idx]
      tst_Z <- Z[tst_idx]
    } else {
      trn_Z <- Z[trn_idx, ]
      tst_Z <- Z[tst_idx, ]
    }

    # From the training set, we fit FPCA decompositions and obtain the parameter estimates #
    Fit <- modelFit_VS(trn_Y, trn_X, trn_Z, startT, pve, nbasis, trn_weight)
    list_fpcaFit <- Fit$list_fpcaFit
    list_cvfit <- Fit$list_cvfit

    for (j in 1:length(tau.sequence)) {
      curr.tau <- tau.sequence[j]
      # use the fitted model to make the prediction for individuals in the test set #
      result <- modelPredict_VS(list_fpcaFit, list_cvfit, tst_X, tst_Z, startT, curr.tau)
      error.matrix[i, j] <- sum(result$final.label != tst_Y, na.rm = TRUE) / sum(!is.na(tst_Y))
      cost.matrix[i, j] <- result$avg.cost
    }
  }

  test.error <- colMeans(error.matrix)
  test.cost <- colMeans(cost.matrix)

  index.cost <- which(test.cost <= budget)                  # index set for tau that satisfies budget constraint

  tau.sequence.satisfied <- tau.sequence[index.cost]
  test.error.satisfied <- test.error[index.cost]
  test.cost.satisfied <- test.cost[index.cost]

  final.tau <- tau.sequence.satisfied[max(which(test.error.satisfied == min(test.error.satisfied)))]

  final.result <- modelFit_VS(Y, X, Z, startT, pve, nbasis, weight)

  return (list("final.result" = final.result, "final.tau" = final.tau))
}

#' Example training data (high dimensional Z)
#'
#' @format A data frame with 400 rows and 112 columns. The 1st column is the outcome variable Y, starting from the 2nd
#' column to 62nd column is the longitudinal biomarker at 61 time grids, the 63rd column to 112nd column are other baseline covariate Z.
#'
#' @source A simulated data set
"train_data_mulZ"

#' Example test data (high dimensional Z)
#'
#' @format A data frame with 2000 rows and 112 columns. The 1st column is the outcome variable Y, starting from the 2nd
#' column to 62nd column is the longitudinal biomarker at 61 time grids, the 63rd column to 112nd column are other baseline covariate Z.
#'
#' @source A simulated data set
"test_data_mulZ"
