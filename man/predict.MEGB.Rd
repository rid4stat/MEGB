% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/MEGB.R
\name{predict.MEGB}
\alias{predict.MEGB}
\title{Predict with longitudinal trees and random forests.}
\usage{
\method{predict}{MEGB}(object, X, Z, id, time, ntree, ...)
}
\arguments{
\item{object}{: a \code{longituRF} output of (S)MERF; (S)REEMforest; (S)MERT or (S)REEMtree function.}

\item{X}{[matrix]: matrix of the fixed effects for the new observations to be predicted.}

\item{Z}{[matrix]: matrix of the random effects for the new observations to be predicted.}

\item{id}{[vector]: vector of the identifiers of the new observations to be predicted.}

\item{time}{[vector]: vector of the time measurements of the new observations to be predicted.}

\item{ntree}{[numeric]: Number of trees to be used in prediction not less than number of trees used in the model object MEGB. The default value is \code{ntree=500}.}

\item{...}{: low levels arguments.}
}
\value{
vector of the predicted output for the new observations.
}
\description{
Predict with longitudinal trees and random forests.
}
\examples{
set.seed(1)
data <-simLong(n = 20,p = 6,rel_p = 6,time_points = 10,rho_W = 0.6, rho_Z=0.6,
              random_sd_intercept = sqrt(0.5),
              random_sd_slope = sqrt(3),
              noise_sd = 0.5,linear=TRUE)  # Generate the data composed by n=20 individuals.
# Train a MEGB model on the generated data. Should take ~ 7 seconds
megb <-   MEGB(X=as.matrix(data[,-1:-5]),Y=as.matrix(data$Y),
Z=as.matrix(data[,4:5]),id=data$id,time=data$time,ntree=500,verbose=TRUE)
# Then we predict on the learning sample :
pred.MEGB <- predict(megb, X=as.matrix(data[,-1:-5]), Z=as.matrix(data[,4:5]),
id=data$id, time=data$time,ntree=500)
# Let's have a look at the predictions
# the predictions are in red while the real output trajectories are in blue:
par(mfrow=c(4,5),mar=c(2,2,2,2))
for (i in unique(data$id)){
  w <- which(data$id==i)
  plot(data$time[w],data$Y[w],type="l",col="blue")
  lines(data$time[w],pred.MEGB[w], col="red")
}

}
