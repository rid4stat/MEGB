# MEGB
MEGB: An R package for Mixed Effect Gradient Boosting for High-dimensional Longitudinal Data
# MEGB: Mixed Effect Gradient Boosting Algorithm

MEGB is an adaptation of the gradient boosting regression method tailored for longitudinal data, inspired by the Mixed Effect Random Forest (MERF) model developed by Hajjem et al. (2014). The algorithm estimates the parameters of a semi-parametric mixed-effects model:

\[ Y_i(t) = f(X_i(t)) + Z_i(t)\beta_i + \epsilon_i \]

where:
- \( Y_i(t) \): Output at time \( t \) for the \( i \)-th individual.
- \( X_i(t) \): Predictors for fixed effects at time \( t \).
- \( Z_i(t) \): Predictors for random effects at time \( t \).
- \( \beta_i \): Random effect coefficients.
- \( \epsilon_i \): Residual error.

## Usage

```R
MEGB(X, Y, id, Z, iter = 100, ntree = 500, time, shrinkage = 0.05,
     interaction.depth = 1, n.minobsinnode = 5, cv.folds = 0, delta = 0.001, verbose = TRUE)

