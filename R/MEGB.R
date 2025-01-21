#' Mixed Effect Gradient Boosting (MEGB) Algorithm
#'
#' MEGB is an adaptation of the gradient boosting regression method to longitudinal data similar to the Mixed Effect Random Forest (MERF) developed by Hajjem et. al. (2014) <doi:10.1080/00949655.2012.741599> which was implemented by Capitaine et. al. (2020) <doi:10.1177/0962280220946080>.
#' The algorithm estimates the parameters of a semi-parametric mixed-effects model: \deqn{Y_i(t)=f(X_i(t))+Z_i(t)\beta_i+\epsilon_i}
#' with \eqn{Y_i(t)} the output at time \eqn{t} for the \eqn{i}th individual; \eqn{X_i(t)} the input predictors (fixed effects) at time \eqn{t} for the \eqn{i}th individual;
#' \eqn{Z_i(t)} are the random effects at time \eqn{t} for the \eqn{i}th individual;
#'  \eqn{\epsilon_i} is the residual error.
#'
#' @param X [matrix]: A \code{N} x \code{p} matrix containing the \code{p} predictors of the fixed effects, column codes for a predictor.
#' @param Y [vector]: A vector containing the output trajectories.
#' @param id [vector]: Is the vector of the identifiers for the different trajectories.
#' @param Z [matrix]: A \code{N} x \code{q} matrix containing the \code{q} predictor of the random effects.
#' @param iter [numeric]: Maximal number of iterations of the algorithm. The default is set to \code{iter=100}
#' @param ntree [numeric]: Number of trees to grow. This should not be set to too small a number, to ensure that every input row gets predicted at least a few times. The default value is \code{ntree=500}.
#' @param time [vector]: Is the vector of the measurement times associated with the trajectories in \code{Y},\code{Z} and \code{X}.
#' @param shrinkage [numeric]: a shrinkage parameter applied to each tree in the expansion. Also known as the learning rate or step-size reduction. The default value is set to 0.05.
#' @param interaction.depth [numeric]: The maximum depth of variable interactions: 1 builds an additive model, 2 builds a model with up to two-way interactions, etc. The default value is set to 1.
#' @param n.minobsinnode [numeric]: minimum number of observations (not total weights) in the terminal nodes of the trees. The default value is set to 5.
#' @param cv.folds [numeric]: Number of cross-validation folds to perform. If cv.folds>1 then gbm, in addition to the usual fit, will perform a cross-validation and calculate an estimate of generalization error returned in cv_error. The default value is set to 0.
#' @param delta [numeric]: The algorithm stops when the difference in log likelihood between two iterations is smaller than \code{delta}. The default value is set to 0.001
#' @param verbose [boolean]: If TRUE, MEGB will print out number of iterations to achieve convergence. Default is TRUE.
#'
#' @import gbm
#' @import stats
#' @return A fitted MEGB model which is a list of the following elements: \itemize{
#' \item \code{forest:} GBMFit obtained at the last iteration.
#' \item \code{random_effects :} Predictions of random effects for different trajectories.
#' \item \code{id_btilde:} Identifiers of individuals associated with the predictions \code{random_effects}.
#' \item \code{var_random_effects: } Estimation of the variance covariance matrix of random effects.
#' \item \code{sigma: } Estimation of the residual variance parameter.
#' \item \code{time: } The vector of the measurement times associated with the trajectories in \code{Y},\code{Z} and \code{X}.
#' \item \code{LL:} Log-likelihood of the different iterations.
#' \item \code{id: } Vector of the identifiers for the different trajectories.
#' \item \code{OOB: } OOB error of the fitted random forest at each iteration.
#' }
#'
#' @export
#'
#' @examples
#' set.seed(1)
#' data <-simLong(n = 20,p = 6,rel_p = 6,time_points = 10,rho_W = 0.6, rho_Z=0.6,
#'               random_sd_intercept = sqrt(0.5),
#'               random_sd_slope = sqrt(3),
#'               noise_sd = 0.5,linear=TRUE)  # Generate the data composed by n=20 individuals.
#' # Train a MEGB model on the generated data. Should take ~ 7 seconds
#' megb <-   MEGB(X=as.matrix(data[,-1:-5]),Y=as.matrix(data$Y),
#' Z=as.matrix(data[,4:5]),id=data$id,time=data$time,ntree=500,cv.folds=0,verbose=TRUE)
#' megb$forest # is the fitted gradient boosting (GBMFit) (obtained at the last iteration).
#' megb$random_effects # are the predicted random effects for each individual.
#' plot(megb$LL,type="o",col=2) # evolution of the log-likelihood.
#' megb$OOB # OOB error at each iteration.
#'
#'
MEGB = function (X, Y, id, Z, iter = 100, ntree = 500, time, shrinkage=0.05,
                 interaction.depth=1, n.minobsinnode=5, cv.folds=0, delta = 0.001,verbose=TRUE)
{


  q <- dim(Z)[2]
  nind <- length(unique(id))
  btilde <- matrix(0, nind, q)
  sigmahat <- 1
  Btilde <- diag(rep(1, q))
  epsilonhat <- rep(0, length(Y))
  id_btilde <- unique(id)
  Tiime <- sort(unique(time))
  sigma2 <- 1
  Vrai <- NULL
  inc <- 1
  OOB <- NULL


  for (i in 1:iter) {
    ystar <- rep(NA, length(Y))
    for (k in 1:nind) {
      indiv <- which(id == unique(id)[k])
      ystar[indiv] <- Y[indiv] - Z[indiv, , drop = FALSE] %*%
        btilde[k, ]
    }
    set.seed(i)
    forest <- gbm(ystar ~ ., data = data.frame(X,ystar),distribution = "gaussian",
                  n.trees = ntree, shrinkage = shrinkage,
                  interaction.depth = interaction.depth, bag.fraction = 1, train.fraction = 1,
                  n.minobsinnode = n.minobsinnode, cv.folds = cv.folds, keep.data = TRUE,
                  verbose = FALSE)

    fhat <- as.numeric(predict(forest,n.trees=ntree,newdata=data.frame(X,ystar)))
    OOB[i] <- forest$train.error[ntree]
    for (k in 1:nind) {
      indiv <- which(id == unique(id)[k])
      V <- Z[indiv, , drop = FALSE] %*% Btilde %*%
        t(Z[indiv, , drop = FALSE]) + diag(as.numeric(sigmahat),
                                           length(indiv), length(indiv))
      btilde[k, ] <- Btilde %*% t(Z[indiv, , drop = FALSE]) %*%
        solve(V) %*% (Y[indiv] - fhat[indiv])
      epsilonhat[indiv] <- Y[indiv] - fhat[indiv] -
        Z[indiv, , drop = FALSE] %*% btilde[k, ]
    }
    sigm <- sigmahat
    sigmahat <- sig(sigma = sigmahat, id = id, Z = Z,
                    epsilon = epsilonhat, Btilde = Btilde)
    Btilde <- bay(bhat = btilde, Bhat = Btilde, Z = Z,
                  id = id, sigmahat = sigm)
    Vrai <- c(Vrai, logV(Y, fhat, Z, time, id, Btilde,
                         0, sigmahat))
    if (i > 1)
      inc <- abs((Vrai[i - 1] - Vrai[i])/Vrai[i -
                                                1])
    if (inc < delta) {
      if(verbose){print(paste0("stopped after ", i, " iterations."))}
      results <- list(forest = forest, random_effects = btilde,
                      var_random_effects = Btilde, sigma = sigmahat,
                      id_btilde = unique(id), LL = Vrai,
                      id = id, time = time, OOB = OOB)
      class(results) <- "MEGB"
      return(results)
    }
  }

}


#' Predict with longitudinal trees and random forests.
#'
#' @param object : a \code{longituRF} output of (S)MERF; (S)REEMforest; (S)MERT or (S)REEMtree function.
#' @param X [matrix]: matrix of the fixed effects for the new observations to be predicted.
#' @param Z [matrix]: matrix of the random effects for the new observations to be predicted.
#' @param id [vector]: vector of the identifiers of the new observations to be predicted.
#' @param time [vector]: vector of the time measurements of the new observations to be predicted.
#' @param ntree [numeric]: Number of trees to be used in prediction not less than number of trees used in the model object MEGB. The default value is \code{ntree=500}.
#' @param ... : low levels arguments.
#'
#' @import stats
#' @import gbm
#'
#' @return vector of the predicted output for the new observations.
#'
#' @export
#'
#' @examples
#' set.seed(1)
#' data <-simLong(n = 20,p = 6,rel_p = 6,time_points = 10,rho_W = 0.6, rho_Z=0.6,
#'               random_sd_intercept = sqrt(0.5),
#'               random_sd_slope = sqrt(3),
#'               noise_sd = 0.5,linear=TRUE)  # Generate the data composed by n=20 individuals.
#' # Train a MEGB model on the generated data. Should take ~ 7 seconds
#' megb <-   MEGB(X=as.matrix(data[,-1:-5]),Y=as.matrix(data$Y),
#' Z=as.matrix(data[,4:5]),id=data$id,time=data$time,ntree=500,verbose=TRUE)
#' # Then we predict on the learning sample :
#' pred.MEGB <- predict(megb, X=as.matrix(data[,-1:-5]), Z=as.matrix(data[,4:5]),
#' id=data$id, time=data$time,ntree=500)
#' # Let's have a look at the predictions
#' # the predictions are in red while the real output trajectories are in blue:
#' par(mfrow=c(4,5),mar=c(2,2,2,2))
#' for (i in unique(data$id)){
#'   w <- which(data$id==i)
#'   plot(data$time[w],data$Y[w],type="l",col="blue")
#'   lines(data$time[w],pred.MEGB[w], col="red")
#' }
#'
predict.MEGB <- function(object,X,Z,id,time,ntree,...){
  dfX = data.frame(X)
  colnames(dfX) = colnames(X)
  n <- length(unique(id))
  id_btilde <- object$id_btilde
  f <- as.numeric(predict(object$forest,newdata=dfX,n.trees=ntree))
  Time <- object$time
  id_btilde <- unique(id)
  Ypred <- rep(0,length(id))
  id.app=object$id
  for (i in 1:length(unique(id))){
    w <- which(id==unique(id)[i])
    k <- which(id_btilde==unique(id)[i])
    Ypred[w] <- f[w] + Z[w,, drop=FALSE]%*%object$random_effects[k,]
  }
  return(Ypred)

}


#' Title
#'
#'
#' @import stats
#'
#' @keywords internal
sig <- function(sigma,id,Z, epsilon, Btilde){ #### fonction d'actualisation du param?tre de la variance des erreurs
  nind <- length(unique(id))
  Nombre <- length(id)
  sigm <- 0
  for (j in 1:nind){
    w <- which(id==unique(id)[j])
    V <- Z[w,, drop=FALSE]%*%Btilde%*%t(Z[w,, drop=FALSE])+diag(as.numeric(sigma),length(w),length(w))
    sigm <- sigm + t(epsilon[w])%*%epsilon[w] + sigma*(length(w)-sigma*(sum(diag(solve(V)))))
  }
  sigm <- sigm/Nombre
  return(sigm)
}

#' Title
#'
#' @import stats
#'
#' @keywords internal
bay <- function(bhat,Bhat,Z,id, sigmahat){ #### actualisation des param?tres de B
  nind <- length(unique(id))
  q <- dim(Z)[2]
  Nombre <- length(id)
  D <- 0
  for (j in 1:nind){
    w <- which(id==unique(id)[j])
    V <- Z[w,, drop=FALSE]%*%Bhat%*%t(Z[w,, drop=FALSE])+diag(as.numeric(sigmahat),length(w),length(w))
    D <- D+ (bhat[j,]%*%t(bhat[j,]))+ (Bhat- Bhat%*%t(Z[w,, drop=FALSE])%*%solve(V)%*%Z[w,, drop=FALSE]%*%Bhat)
  }
  D <- D/nind
  return(D)
}

#' Title
#'
#' @import stats
#'
#' @keywords internal
logV <- function(Y,f,Z,time,id,B,gamma,sigma){ # Maximization of variance of Y at M stage of EM
  Vraisem <- 0
  for (i in 1:length(unique(id))){
    w <- which(id==unique(id)[i])
    V <- Z[w,,drop=FALSE]%*%B%*%t(Z[w,,drop=FALSE])+diag(as.numeric(sigma),length(w),length(w))
    Vraisem <- Vraisem + log(det(V))+ t(Y[w]-f[w])%*%solve(V)%*%(Y[w]-f[w])
  }
  return(Vraisem)
}




#' Simulate Low/High Dimensional and Linear/Nonlinear Longitudinal dataset.
#'
#'
#' Simulate p-dimensional linear/Nonlinear mixed-effects model given by: \deqn{Y_i(t)=f(X_i(t))+Z_i(t)\beta_i+\epsilon_i}
#' with \eqn{Y_i(t)} the output at time \eqn{t} for the \eqn{i}th individual; \eqn{X_i(t)} the input predictors (fixed effects) at time \eqn{t} for the \eqn{i}th individual;
#' \eqn{Z_i(t)} are the random effects at time \eqn{t} for the \eqn{i}th individual; \eqn{\epsilon_i} is the residual error with
#' variance \eqn{\sigma^2}. If linear, \eqn{f(X_i(t)) = X_i(t)\theta}, where \eqn{\theta = 1, \forall p}, otherwise if nonlinear, the
#' approach by Capitaine et al. (2021) is adapted.
#'
#' @param n [numeric]: Number of individuals.
#' @param p [numeric]: Number of predictors.
#' @param rel_p [numeric]: Number of relevant predictors (true predictors that are correlated to the outcome.). The default value is \code{rel_p=6} if linear and \code{rel_p=2} if nonlinear.
#' @param time_points [numeric]: Number of realizations per individual.  The default value is \code{time_points=10}.
#' @param rho_W [numeric]: Within subject correlation. The default value is \code{rho_W=0.5}.
#' @param rho_Z [numeric]: Correlation between intercept and slope for the random effect coefficients. The default value is \code{rho_Z=0.5}.
#' @param random_sd_intercept [numeric]: Standard deviation for the random intercept. The default value is \code{random_sd_intercept=}\eqn{\sqrt{0.5}}.
#' @param random_sd_slope [numeric]: Standard deviation for the random slope. The default value is \code{random_sd_slope=}\eqn{\sqrt{3}}.
#' @param noise_sd [numeric]: Standard deviation for the random slope. The default value is \code{noise_sd=0.5}.
#' @param linear [boolean]: If TRUE, a linear mixed effect model is simulated, if otherwise, a semi-parametric model similar to the one used in Capitaine et al. (2021).
#'
#' @import MASS
#' @import latex2exp
#'
#' @return a dataframe of dimension (n*time_points) by (p+5) containing the following elements: \itemize{
#' \item \code{id:} vector of the individual IDs.
#' \item \code{time:} vector of the time realizations.
#' \item \code{Y:} vector of the outcomes variable.
#' \item \code{RandomIntercept:} vector of the Random Intercept.
#' \item \code{RandomSlope:} vector of the Random Slope.
#' \item \code{Vars :} Remainder columns corresponding to the fixed effect variables.
#' }
#'
#' @export
#'
#' @examples
#' set.seed(1)
#' data = simLong(n = 17,p = 6,rel_p = 6,time_points = 10,rho_W = 0.6, rho_Z=0.6,
#'               random_sd_intercept = sqrt(0.5),
#'               random_sd_slope = sqrt(3),
#'               noise_sd = 0.5,linear=FALSE) # Generate the data
#' head(data)   # first six rows of the data.
#' # Let's see the output :
#' w <- which(data$id==1)
#' plot(data$time[w],data$Y[w],type="l",ylim=c(min(data$Y),max(data$Y)), col="grey")
#' for (i in unique(data$id)){
#'   w <- which(data$id==i)
#'   lines(data$time[w],data$Y[w], col='grey')
#' }
#' # Let's see the fixed effects predictors:
#' par(mfrow=c(2,3), mar=c(2,3,3,2))
#' for (i in 1:ncol(data[,-1:-5])){
#'   w <- which(data$id==1)
#'   plot(data$time[w],data[,-1:-5][w,i], col="grey",ylim=c(min(data[,-1:-5][,i]),
#'   max(data[,-1:-5][,i])),xlim=c(1,max(data$time)),main=latex2exp::TeX(paste0("$X^{(",i,")}$")))
#'   for (k in unique(data$id)){
#'     w <- which(data$id==k)
#'     lines(data$time[w],data[,-1:-5][w,i], col="grey")
#'   }
#' }
#'
simLong <- function(
    n, p, rel_p = 6, time_points, rho_W = 0.5, rho_Z = 0.5,
    random_sd_intercept = 2, random_sd_slope = 1, noise_sd = 1,linear=T
) {

  time <- rep(1:time_points,n)
  id = rep(1:n,each=time_points)


  # Generate covariance matrix with AR(1) structure for within-subject correlation
  Sigma_within <- matrix(0, nrow = time_points, ncol = time_points)
  for (i in 1:time_points) {
    for (j in 1:time_points) {
      Sigma_within[i, j] <- rho_W^abs(i - j)
    }
  }

  # Fixed effect predictors simulation

  Xf = list()

  for (i in 1:n) {
    # Generate predictors with AR(1) correlation structure
    Xp = matrix(NA,nrow=time_points,ncol=p)
    for (j in 1:p) {
      Xp[, j] <- MASS::mvrnorm(1, mu = rep(0, time_points), Sigma = Sigma_within)
    }
    Xf[[i]] = Xp
  }

  X = t(matrix(sapply(Xf,function(x)x),p,n*time_points))


  # Generate random intercepts and slopes for subjects

  SigmaZ = matrix(c(random_sd_intercept^2,rho_Z,rho_Z,random_sd_slope^2),2,2)

  b<- matrix(0,length(unique(id)),2)
  for (i in 1:length(unique(id))){
    b[i,] <- MASS::mvrnorm(1, mu = rep(0, 2), Sigma = SigmaZ)
  }

  Z <- as.matrix(cbind(rep(1,length(id)),2*runif(length(id))))

  random_effects  <- NULL
  for (i in 1:length(unique(id))){
    w <- which(id==unique(id)[i])
    random_effects <- c(random_effects, Z[w,, drop=FALSE]%*%b[i,])
  }

  # Add response variable with mixed-effects model structure
  beta <- rep(1,rel_p)  # Fixed effects coefficients if Linear
  noise <- rnorm(n, mean = 0, sd = noise_sd)  # Random noise

  if(linear){

    Y <- X[,1:rel_p] %*% beta + random_effects + noise  # Response variable
  }else{

    # Below Nonlinear fixed effect predictors simulation was adapted from Capitaine et al. (2021)

    X[,1] <- 2.44+0.04*(time-((time-6)^2)/(time/3))+rnorm(n,0,0.2)
    X[,2] <- 0.5*time-0.1*(time-5)^2+rnorm(n,0,0.2)
    X[,3] <- 0.25*time-0.05*(time-6)^2+rnorm(n,0,0.2)
    X[,4] <- cos((time-1)/3)+rnorm(n,0,0.2)
    X[,5] <- 0.1*time + sin(0.6*time+1.3)+rnorm(n,0,0.2)
    X[,6] <- -0.1*time^2+rnorm(n,0,0.2)

    #Friedman 5-dimensional
    #fixed_effect = 10*sin(pi*X[,1]*X[,2])+20*(X[,3]-0.5)^2+10*X[,4]+5*X[,5]+X[,6]

    fixed_effect = 1.3*(X[,1])^2+2*abs(X[,2])^0.5 # +2*(X[,3]-0.5)^2 + 2*abs(X[,4])+X[,5]+X[,6]

    Y <- fixed_effect + random_effects + noise  # Response variable
  }


  # Combine data, including random effects used
  temp_data <- data.frame(
    id = id,
    time = time,
    Y = as.vector(Y),
    RandomIntercept = Z[,1],
    RandomSlope = Z[,2],
    as.data.frame(X)  # Predictors
  )
  data <- temp_data


  # Rename predictor columns
  colnames(data)[5 + (1:p)] <- paste0("Var", 1:p)

  return(data)
}

