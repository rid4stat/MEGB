% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/MEGB.R
\name{MEGB}
\alias{MEGB}
\title{Mixed Effect Gradient Boosting (MEGB) Algorithm}
\usage{
MEGB(
  X,
  Y,
  id,
  Z,
  iter = 100,
  ntree = 500,
  time,
  shrinkage = 0.05,
  interaction.depth = 1,
  n.minobsinnode = 5,
  cv.folds = 0,
  delta = 0.001,
  verbose = TRUE
)
}
\arguments{
\item{X}{[matrix]: A \code{N} x \code{p} matrix containing the \code{p} predictors of the fixed effects, column codes for a predictor.}

\item{Y}{[vector]: A vector containing the output trajectories.}

\item{id}{[vector]: Is the vector of the identifiers for the different trajectories.}

\item{Z}{[matrix]: A \code{N} x \code{q} matrix containing the \code{q} predictor of the random effects.}

\item{iter}{[numeric]: Maximal number of iterations of the algorithm. The default is set to \code{iter=100}}

\item{ntree}{[numeric]: Number of trees to grow. This should not be set to too small a number, to ensure that every input row gets predicted at least a few times. The default value is \code{ntree=500}.}

\item{time}{[vector]: Is the vector of the measurement times associated with the trajectories in \code{Y},\code{Z} and \code{X}.}

\item{shrinkage}{[numeric]: a shrinkage parameter applied to each tree in the expansion. Also known as the learning rate or step-size reduction. The default value is set to 0.05.}

\item{interaction.depth}{[numeric]: The maximum depth of variable interactions: 1 builds an additive model, 2 builds a model with up to two-way interactions, etc. The default value is set to 1.}

\item{n.minobsinnode}{[numeric]: minimum number of observations (not total weights) in the terminal nodes of the trees. The default value is set to 5.}

\item{cv.folds}{[numeric]: Number of cross-validation folds to perform. If cv.folds>1 then gbm, in addition to the usual fit, will perform a cross-validation and calculate an estimate of generalization error returned in cv_error. The default value is set to 0.}

\item{delta}{[numeric]: The algorithm stops when the difference in log likelihood between two iterations is smaller than \code{delta}. The default value is set to 0.001}

\item{verbose}{[boolean]: If TRUE, MEGB will print out number of iterations to achieve convergence. Default is TRUE.}
}
\value{
A fitted MEGB model which is a list of the following elements: \itemize{
\item \code{forest:} GBMFit obtained at the last iteration.
\item \code{random_effects :} Predictions of random effects for different trajectories.
\item \code{id_btilde:} Identifiers of individuals associated with the predictions \code{random_effects}.
\item \code{var_random_effects: } Estimation of the variance covariance matrix of random effects.
\item \code{sigma: } Estimation of the residual variance parameter.
\item \code{time: } The vector of the measurement times associated with the trajectories in \code{Y},\code{Z} and \code{X}.
\item \code{LL:} Log-likelihood of the different iterations.
\item \code{id: } Vector of the identifiers for the different trajectories.
\item \code{OOB: } OOB error of the fitted random forest at each iteration.
}
}
\description{
MEGB is an adaptation of the gradient boosting regression method to longitudinal data similar to the Mixed Effect Random Forest (MERF) developed by Hajjem et. al. (2014) <doi:10.1080/00949655.2012.741599> which was implemented by Capitaine et. al. (2020) <doi:10.1177/0962280220946080>.
The algorithm estimates the parameters of a semi-parametric mixed-effects model: \deqn{Y_i(t)=f(X_i(t))+Z_i(t)\beta_i+\epsilon_i}
with \eqn{Y_i(t)} the output at time \eqn{t} for the \eqn{i}th individual; \eqn{X_i(t)} the input predictors (fixed effects) at time \eqn{t} for the \eqn{i}th individual;
\eqn{Z_i(t)} are the random effects at time \eqn{t} for the \eqn{i}th individual;
 \eqn{\epsilon_i} is the residual error.
}
\examples{
set.seed(1)
data <-simLong(n = 20,p = 6,rel_p = 6,time_points = 10,rho_W = 0.6, rho_Z=0.6,
              random_sd_intercept = sqrt(0.5),
              random_sd_slope = sqrt(3),
              noise_sd = 0.5,linear=TRUE)  # Generate the data composed by n=20 individuals.
# Train a MEGB model on the generated data. Should take ~ 7 seconds
megb <-   MEGB(X=as.matrix(data[,-1:-5]),Y=as.matrix(data$Y),
Z=as.matrix(data[,4:5]),id=data$id,time=data$time,ntree=500,cv.folds=0,verbose=TRUE)
megb$forest # is the fitted gradient boosting (GBMFit) (obtained at the last iteration).
megb$random_effects # are the predicted random effects for each individual.
plot(megb$LL,type="o",col=2) # evolution of the log-likelihood.
megb$OOB # OOB error at each iteration.


}
