% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/reinforced_VS.R
\name{modelFit_VS}
\alias{modelFit_VS}
\title{Model fit for the training set, variable selection version}
\usage{
modelFit_VS(Y, X, Z, startT, pve, nbasis, weight)
}
\arguments{
\item{Y}{The outcome variable, vector of length \eqn{n}, taking values in \eqn{{1, 0, NA}}, where 1 = disease, 0 = not, NA = missing.}

\item{X}{Observed longitudinal biomarker, matrix of \eqn{n} by \eqn{nTotal}, where \eqn{nTotal} denotes the total number of time grids.
Missing values are denoted by NA.}

\item{Z}{Other baseline covariates.}

\item{startT}{Time of the first prediction, denoted by \eqn{t_1} in the manuscript. For instance, if the time grids are \eqn{{0,1/60,2/60,...,1}},
then startT = 25 means that the first prediction is made at \eqn{t = 24/60}.}

\item{pve}{Proportion of variance explained in FPCA.}

\item{nbasis}{Number of B-spline basis functions needed for estimation of the mean function and smoothing of covariance.}

\item{weight}{Weight for each individual.}
}
\value{
\item{list_fpcaFit}{FPCA decomposition at each time grid from startT to the end.}
\item{list_cvfit}{Elastic net logistic regression at each time grid from startT to the end.}
}
\description{
\code{modelFit_VS} outputs the FPCA (functional principal component analysis) decomposition and the elastic net logistic regression
at each time grid. This function is used when the baseline covariates Z are high-dimensional.
}
\examples{
\donttest{
library(reinforcedPred)

# take the example training data (high dimensional Z) from the reinforcedPred package
# see documentation for details about the data set train_data_mulZ
Y <- as.numeric(train_data_mulZ$Y)
tildeX.missing <- as.matrix(train_data_mulZ[,2:62])
Z <- as.matrix(train_data_mulZ[,63:dim(train_data_mulZ)[2]])

# analysis starts
startT <- 25
weight <- rep(1, length(Y))

result <- modelFit_VS(Y, tildeX.missing, Z, startT, pve = 0.99, nbasis = 10, weight)

# obtained elastic net logistic regression fit and FPCA decompositions
list_cvfit <- result$list_cvfit
list_fpcaFit <- result$list_fpcaFit
}
}
