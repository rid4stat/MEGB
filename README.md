# MEGB
MEGB: An R package for Mixed Effect Gradient Boosting for High-dimensional Longitudinal Data
# Mixed Effect Gradient Boosting (MEGB) Algorithm

MEGB is an adaptation of the gradient boosting regression method to longitudinal data, similar to the Mixed Effect Random Forest (MERF) developed by Hajjem et al. (2014) and implemented by Capitaine et al. (2020). The algorithm estimates the parameters of a semi-parametric mixed-effects model:

$$ Y_i(t) = f(X_i(t)) + Z_i(t)\beta_i + \epsilon_i $$

where:
- $Y_i(t)$: Output at time $t$ for the $i$-th individual.
- $X_i(t)$: Input predictors (fixed effects) at time $t$ for the $i$-th individual.
- $Z_i(t)$: Random effects at time $t$ for the $i$-th individual.
- $\epsilon_i$: Residual error.

## Usage

```R
MEGB(
  X,       # [matrix] Predictors of fixed effects (N x p matrix).
  Y,       # [vector] Output trajectories.
  id,      # [vector] Identifiers for trajectories.
  Z,       # [matrix] Predictors of random effects (N x q matrix).
  iter = 100,             # [numeric] Maximum number of iterations (default: 100).
  ntree = 500,            # [numeric] Number of trees to grow (default: 500).
  time,                   # [vector] Measurement times for trajectories.
  shrinkage = 0.05,       # [numeric] Learning rate (default: 0.05).
  interaction.depth = 1,  # [numeric] Maximum depth of variable interactions (default: 1).
  n.minobsinnode = 5,     # [numeric] Minimum observations in terminal nodes (default: 5).
  cv.folds = 0,           # [numeric] Number of cross-validation folds (default: 0).
  delta = 0.001,          # [numeric] Convergence threshold (default: 0.001).
  verbose = TRUE          # [boolean] Verbose output (default: TRUE).
)
```

## Output

A fitted MEGB model, which includes:

- `forest`: GBMFit obtained at the last iteration.
- `random_effects`: Predictions of random effects for trajectories.
- `id_btilde`: Identifiers of individuals associated with `random_effects`.
- `var_random_effects`: Variance-covariance matrix of random effects.
- `sigma`: Residual variance parameter estimate.
- `time`: Measurement times for trajectories.
- `LL`: Log-likelihood across iterations.
- `id`: Identifiers for trajectories.
- `OOB`: Out-of-bag error for the fitted gradient boosting at each iteration.

## Examples

### Training an MEGB Model

```R
# MEGB can be installed from CRAN:
install.packages("MEGB")
# The development version can be installed from GitHub:
devtools::install_github("rid4stat/MEGB")
library(MEGB)
set.seed(1)
data <- simLong(
  n = 20, p = 6, rel_p = 6, time_points = 10, rho_W = 0.6, rho_Z = 0.6,
  random_sd_intercept = sqrt(0.5), random_sd_slope = sqrt(3),
  noise_sd = 0.5, linear = TRUE
)  # Generate data with n = 20 individuals.

megb <- MEGB(
  X = as.matrix(data[,-1:-5]),
  Y = as.matrix(data$Y),
  Z = as.matrix(data[,4:5]),
  id = data$id,
  time = data$time,
  ntree = 500,
  cv.folds = 0,
  verbose = TRUE
)

# View results
megb$forest        # Fitted gradient boosting (GBMFit).
megb$random_effects # Predicted random effects for each individual.
plot(megb$LL, type = "o", col = 2) # Evolution of log-likelihood.
megb$OOB           # Out-of-bag error at each iteration.
```

### Making Predictions

```R
pred.MEGB <- predict(
  megb,
  X = as.matrix(data[,-1:-5]),
  Z = as.matrix(data[,4:5]),
  id = data$id,
  time = data$time,
  ntree = 500
)

# Visualizing predictions
par(mfrow = c(4, 5), mar = c(2, 2, 2, 2))
for (i in unique(data$id)) {
  w <- which(data$id == i)
  plot(data$time[w], data$Y[w], type = "l", col = "blue")
  lines(data$time[w], pred.MEGB[w], col = "red")
}
```
# References
- Hajjem, A., Bellavance, F., & Larocque, D. (2014). Mixed-effects random forest for clustered data. Journal of Statistical Computation and Simulation, 84(6), 1313-1328.
- Capitaine, L., Genuer, R., & Thiébaut, R. (2021). Random forests for high-dimensional longitudinal data. Statistical methods in medical research, 30(1), 166-184.
