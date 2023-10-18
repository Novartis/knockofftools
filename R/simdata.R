#' Simulate Gaussian and binary covariate predictors
#'
#' simulate Gaussian predictors with mean zero
#' and covariance structure determined by "cov_type" argument. Then p_b
#' randomly selected columns are dichotomized.
#'
#' @param n number of observations (rows of X)
#' @param p total number of covariates (columns of X) both continuous and binary
#' @param p_b number of binary covariates (0 <= p_b <= p)
#' @param cov_type character string specifying the covariance function. Can be one of
#' "cov_diag" (independent columns), "cov_equi" (equi-correlated columns), or "cov_ar1" (ar1-correlated columns).
#' The columns are shuffled during simulation
#' @param rho correlation parameter; input to the cov_type function
#'
#' @details This function simulates a data frame, whose rows are multivariate Gaussian with mean zero
#' and covariance structure determined by "cov_type" argument. Then p_b randomly selected columns are
#' dichotomized with the function 1(x>0). The continuous columns are of class "numeric" and the binary
#' columns are set to class "factor".
#'
#' @return the simulated data.frame with n rows and p columns (p_b of which are binary and p-p_b of which are gaussian).
#' Each column is either of class "numeric" or "factor".
#' @export
#'
#' @examples
#' library(knockofftools)
#'
#' # all columns are continuous:
#' X <- generate_X(n=100, p=6, p_b=0, cov_type="cov_equi", rho=0.5)
#'
#' round(cor(X), 2)
#'
#' # two of the six columns are dichotomized (and set to class factor):
#' X <- generate_X(n=100, p=6, p_b=2, cov_type="cov_equi", rho=0.5)
#'
#' # The class of each column:
#' unlist(lapply(X, class))
generate_X <- function(n, p, p_b, cov_type, rho=0.5) {

  sigma_z <- do.call(get(cov_type), args=list(n=n, p=p, rho=rho))

  X <- data.frame(matrix(rnorm(n*p), nrow=n) %*% chol(sigma_z))

  # Random indices for continuous and binary columns of X
  inds_b <- sample(1:p, size=p_b, replace=FALSE)

  # set p_b randomly selected indices as binary factor \delta(Z > 0)
  X[, inds_b] <- lapply((1+sign(X[, inds_b]))/2, as.factor)

  X <- dplyr::mutate_if(X, is.numeric, function(x) as.numeric(scale(x)))

  return(X)

}

# Covariance matrices scaled to be approximately 1/n on the diagonal

cov_diag <- function(n, p, rho=NULL) {
  # Diagonal covariance
  s <- diag(p)
  return(s)
}

cov_equi <- function(n, p, rho = 0.5) {
  # Equicorrelated covariance
  s <- (diag(1 - rho, p, p) + rho)
  return(s)
}

cov_ar1 <- function(n, p, rho = 0.5) {
  # AR(1) covariance
  s <- toeplitz(rho^(0:(p - 1)))
  return(s)
}


#' Simulate Gaussian response from a sparse regression model
#'
#' @param X data.frame with numeric and factor columns only.
#' @param p_nn number of non-null covariate predictors.
#' The regression coefficients (beta) corresponding to
#' columns 1:p_nn of x will be non-zero, all other are set to zero.
#' @param a amplitude of non-null regression coefficients
#'
#' @details This function takes as input data.frame X (created with the function \code{generate_X})
#' that may consist of both numeric and binary factor columns. This data frame is then expanded
#' to a model matrix x (with the model.matrix function) and subsequently scaled in the same way as
#' LASSO scaling. Next we simulate y ~ N(x%*%beta,I) where the first p_nn beta coefficients are equal to a, while
#' the remaining coefficients (p_nn+1):ncol(x) are set to zero.
#'
#' @return simulated Gaussian response from regression model y = x%*%beta+epsilon, where epsilon ~ N(0,I) and
#' x is the (scaled) model.matrix of X.
#'
#' @export
#'
#' @examples
#' library(knockofftools)
#'
#' set.seed(1)
#'
#' # Simulate 4 Gaussian and 2 binary covariate predictors:
#' X <- generate_X(n=100, p=6, p_b=2, cov_type="cov_equi", rho=0.5)
#'
#' # Simulate response from model y = 2*X[,1] + 2*X[,2] + epsilon, where epsilon ~ N(0,1)
#' y <- generate_y(X, p_nn=2, a=2)
generate_y <- function(X, p_nn, a) {

  x <- model.matrix(~., data=X)[,-1]

  n <- nrow(x)
  p <- ncol(x)

  # Standardize design matrix to zero mean, "variance" one
  x_centered <- apply(x, 2, function(x) x - mean(x))
  x <- apply(x_centered, 2, function(x) x / sqrt(mean(x^2)))

  beta <- rep(c(a, 0), c(p_nn, p - p_nn))

  mu <- as.numeric(x %*% beta)

  y <- mu + rnorm(n)

  return(y)

}


#' Function that simulates response from Cox model with Weibull baseline hazard:
#'
#' @param N sample size
#' @param lambda0 baseline hazard scale parameter
#' @param rho baseline hazard shape parameter
#' @param lp linear predictor
#' @export
#'
#' @return survival object with simulated event times under mild censoring
#' @examples
#' # Simulate 10 Gaussian covariate predictors:
#' X <- generate_X(n=100, p=10, p_b=0, cov_type="cov_equi", rho=0.2)
#'
#' # create linear predictor with first 5 beta-coefficients = 1 (all other zero)
#' lp <- generate_lp(X, p_nn = 5, a=1)
#'
#' # Simulate from Weibull hazard with with baseline hazard h0(t) = lambda*rho*t^(rho-1)
#' # and linear predictor, whose first 3 coefficients are non-zero:
#' y <- simulWeib(N=nrow(X), lambda0=0.01, rho=1, lp=lp)
simulWeib <- function(N, lambda0, rho, lp) {

  # Censoring times ~ Exponential(lambdaC)
  lambdaC=0.0005 # very mild censoring
  # Simulate Weibull latent event times
  v <- runif(n=N)
  Tlat <- (- log(v) / (lambda0 * exp(lp)))^(1 / rho)
  # Simulate censoring times
  C <- rexp(n=N, rate=lambdaC)

  # Calculate follow-up times and event indicators
  time <- pmin(Tlat, C)
  status <- as.numeric(Tlat <= C)

  survival::Surv(time=time, event=status)

}


#' Generate linear predictor with first p_nn beta coefficients = a, all other = 0
#'
#' @param X data.frame with numeric and factor columns only
#' @param p_nn number of non-null covariate predictors.
#' The regression coefficients (beta) corresponding to
#' columns 1:p_nn of x will be non-zero, all other are set to zero.
#' @param a amplitude of non-null regression coefficients
#'
#' @return linear predictor X%*%beta
#' @export
#'
#' @keywords internal
generate_lp <- function(X, p_nn, a) {

  x <- model.matrix(~., data=X)[,-1]

  n <- nrow(x)
  p <- ncol(x)

  # Standardize design matrix to zero mean, "variance" one
  x_centered <- apply(x, 2, function(x) x - mean(x))
  x <- apply(x_centered, 2, function(x) x / sqrt(mean(x^2)))

  beta <- rep(c(a, 0), c(p_nn, p - p_nn))

  lp <- as.numeric(x %*% beta)

  return(lp)

}


#' Generate simulated data set for vignettes
#'
#' @return generates the simulated data as obtained by the call data(simdata)
#' @export
#'
#' @examples
#' library(knockofftools)
#' data(simdata)
#' all.equal(generate_simdata(), simdata)
generate_simdata <- function() {

  RNGkind("L'Ecuyer-CMRG")
  set.seed(56969)

  N <- 2000
  p <- 30
  p_b = 10
  p_nn <- 10

  # Generate a 2000 x 30 Gaussian data.frame under equi-correlation(rho=0.5) structure,
  # with 10 of the columns dichotomized
  X <- generate_X(n=N, p=p, p_b=p_b, cov_type = "cov_equi", rho=0.5)

  # Generate linear predictor lp = X%*%beta where first 10 beta-coefficients are = a, all other = 0.
  lp <- generate_lp(X=X, p_nn=p_nn, a=1)

  # simulate gaussian response with mean = lp and sd = 1.
  Yg <- lp + rnorm(N)

  # simulate bernoulli response with mean = exp(lp)/(1+exp(lp)).
  Yb <- factor(rbinom(N, size=1, prob=exp(lp)/(1+exp(lp))))

  # simulate censored survival times from Cox regression with linear predictor lp:
  Tc <- simulWeib(N=N, lambda = 0.01, rho = 1, lp = lp)

  dat <- data.frame(Yg, Yb, Tc, X)

  return(dat)

}
