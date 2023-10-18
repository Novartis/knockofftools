#' Sequential knockoffs for continuous and categorical variables
#'
#' @param X data.frame (or tibble) with "numeric" and "factor" columns only. The number of columns, ncol(X) needs to be > 2.
#' @param seq_simulator function that simulates sequential knockoffs. Default is the function \code{sim_EN}, which simulates response from an estimated elastic-net model
#' @param ... other parameters passed to the function seq_simulator. For the default (elastic-net sequential seq_simulator, \code{seq_simulator = sim_EN})
#' these other parameters are passed to cv.glmnet.
#'
#' @details \code{knockoffs_seq} performs sequential knockoff simulation using elastic-net regression.
#' @return sequential knockoff copy of X. A data.frame or tibble of same type and dimensions as X.
#' @export
#'
#' @examples
#' library(knockofftools)
#'
#' set.seed(1)
#'
#' X <- generate_X(n=100, p=6, p_b=2, cov_type="cov_equi", rho=0.5)
#'
#' # knockoffs based on sequential elastic-net regression:
#' Xk <- knockoffs_seq(X)
knockoffs_seq <- function(X, seq_simulator = sim_EN, ...) {

  check_design(X)

  knockoffs <- X

  # add ".tilde" to column names:
  names(knockoffs) <- paste0(names(knockoffs),".tilde")

  # Randomly shuffle column indices of X:
  shf <- sample(ncol(X))

  # Loop through the columns of input data (in random order)
  loop.count <- 1
  for (i in shf) {

    y <- X[[i]] # i-th column serves as response
    Xp <- X[,-i] # columns[-i] serve as predictors

    if (loop.count > 1) Xp <- cbind(knockoffs[,shf[1:(loop.count-1)]], Xp)

    knockoffs[[i]] <- seq_simulator(y = y, X = Xp, ...)

    loop.count <- loop.count + 1

  }

  # remove ".tilde" from column names:
  names(knockoffs) <- gsub(".tilde","", names(knockoffs))

  return(knockoffs)

}


#' Simulate from elastic-net regression model
#'
#' @param y response vector (either "numeric" or "factor") that gets passed to cv.glmnet
#' @param X data.frame of covariates that are passed to cv.glmnet
#' @param ... other parameters passed to the function cv.glmnet
#'
#' @return simulated response
#' @export
#'
#' @examples
#' library(knockofftools)
#'
#' set.seed(1)
#'
#' X = data.frame(matrix(rnorm(100 * 20), 100, 20))
#' y = X[,1] + rnorm(100)
#'
#' # simulate from elastic-net regression:
#' ysim = sim_EN(y=y, X=X)
#'
#' # simulated versus input response:
#' plot(y, ysim)
sim_EN <- function(y, X, ...) {

  x <- model.matrix(~., data = X)[,-1]

  if (is.factor(y)) {

    classes <- levels(y)

    K <- length(classes)

    gm.cv <- glmnet::cv.glmnet(y=y, x=x, family="multinomial", intercept=TRUE, alpha=0.5, ...)

    # Beta coefficients (excluding intercept)
    beta.coefs <- as.numeric(coef(gm.cv, s = "lambda.min")[[2]])[-1]

    mu <- predict(gm.cv, newx=x, type="response", s="lambda.min")

    mat.multinom <- apply(mu, 1, function(prob) rmultinom(n=1, size=1, prob=prob))

    y.sim <- classes[apply((1:K)*mat.multinom, 2, max)]

    y.sim <- factor(y.sim, levels=classes)

    rmse <- NULL

  } else {

    if(!is.numeric(y)) stop("class(y) needs to be either 'numeric' or 'factor'")

    gm.cv <- glmnet::cv.glmnet(y=y, x=x, family="gaussian", intercept=TRUE, alpha=0.5, ...)

    # Beta coefficients (excluding intercept)
    beta.coefs <- as.numeric(coef(gm.cv, s = "lambda.min"))[-1]

    # columns of predictor matrix corresponding to non-zero beta.coefs:
    non.zero.cols <- which(beta.coefs != 0)

    # Total number of non-zero parameters (including intercept, hence + 1)
    s.lambda = length(non.zero.cols) + 1

    mu <- predict(gm.cv, newx=x, type="response", s="lambda.min")

    rmse = sqrt(sum((y-mu)^2)/(length(y) - s.lambda))

    y.sim <- rnorm(n=length(y), mean=mu, sd=rmse)

  }

  return(y.sim)

}



#' Sparse sequential knockoff generation algorithm
#'
#' This function takes as input a data frame X and returns its sparse sequential knockoff copy.
#' Sparse sequential knockoffs first calculates the adjacency matrix of X (i.e. identifies the
#' zeros/non-zeros of the precision matrix of X). Then it proceeds with the usual
#' sequential knockoffs algorithm, except now each sequential regression only includes
#' covariates that correspond to non-zero elements of the precision matrix of X. This
#' reduces the number of covariates per regression in the original sequential knockoff algorithm.
#' To gain additional speed-up (as compared to sequential knockoffs) we apply least squares
#' as the default method for estimation in each regression, unless the number of
#' covariates exceeds half the number of observations, i.e. p > n/2. In this case
#' we apply the elastic net regularized regression.
#'
#' @param X data.frame (or tibble) with "numeric" and "factor" columns only. The number of columns, ncol(X) needs to be > 2.
#' @param sparsity_estimator name of function that estimates sparsity patter. Default is \code{glasso_adjacency_matrix}.
#' @param seq_simulator name of function that used to estimate the conditional distributions in the sequential steps. Default is the function \code{sim_simple}, which is a least squares fit (continuous variables) or multinomial logistic regression (factor variables) respectively.
#'
#' @return sparse sequential knockoff copy of X. A data.frame or tibble of same type and dimensions as X.
#' @export
#'
#' @examples
#' library(knockofftools)
#'
#' set.seed(1)
#'
#' X <- generate_X(n=100, p=6, p_b=2, cov_type="cov_equi", rho=0.5)
#'
#' # knockoffs based on sequential elastic-net regression:
#' Xk <- knockoffs_sparse_seq(X)
knockoffs_sparse_seq <- function(X, sparsity_estimator=glasso_adjacency_matrix){

  check_design(X)

  knockoffs <- X

  # add ".tilde" to column names:
  names(knockoffs) <- paste0(names(knockoffs),".tilde")

  ## calculate the adjacency matrix
  ## This will determine which (original and knockoff) variables to include as predictors in the sequential regressions
  adjacency.matrix <- sparsity_estimator(X)
  diag(adjacency.matrix) <- 0 # for convenience, assume each X_i is not adjacent with itself

  # Randomly shuffle column indices of X:
  shf <- sample(ncol(X))

  # Loop through the columns of input data (in random order)
  loop.count <- 1
  for (i in shf) {
    y <- X[[i]]
    j <- which(adjacency.matrix[,i] != 0)
    Xp <- X[,j, drop = FALSE]

    if (loop.count > 1) Xp <- cbind(knockoffs[,intersect(shf[1:(loop.count-1)], j)], Xp)

    ## generate knockoffs of this feature
    # (if p < n/2, then use least squares, otherwise use elastic net)
    if (ncol(Xp) < length(y)/2) {
      knockoffs[[i]] <- sim_simple(y = y, X = Xp)
    } else {
      knockoffs[[i]] <- sim_EN(y = y, X = Xp)
    }

    loop.count <- loop.count + 1

  }

  # remove ".tilde" from column names:
  names(knockoffs) <- gsub(".tilde","", names(knockoffs))

  return(knockoffs)

}


#' Simple knockoff generator based on least squares fit (continuous variables) or multinomial logistic regression (factor variables) respectively.
#' If X is empty, knockoffs are sampled from the marginal distribution of y
#'
#' @param y response vector (either "numeric" or "factor") that gets passed to cv.glmnet
#' @param X data.frame of covariates that are passed to cv.glmnet
#'
#' @return simulated response
#' @export
#'
sim_simple <- function(y, X) {

  dataset = cbind(y, X)

  if(is.factor(y)) {

    classes <- levels(y)

    if (length(classes)<2) stop("A categorical feature with less than two different levels was provided. This feature is uninformative and should thus be removed.")

    fit <- nnet::multinom(y ~., data=dataset)
    mu <- predict(fit, type="probs")

    if (length(classes)==2){
      ## If y is binary (i.e. has only 2 classes)
      ## mu only contains the probabilities wrt the reference class,
      ## but this will throw an error when passed to predict.
      ## Thus we complement these probabilities.
      mu <- matrix(c(mu, 1-mu), ncol=2)
    }

    mat.multinom <- apply(mu, 1, function(prob) rmultinom(n=1, size=1, prob=prob))

    y.sim <- factor(classes[apply((1:length(classes))*mat.multinom, 2, max)], levels=classes)

  } else { # if y is numeric

    fit <- lm(y ~ ., data=dataset)
    mu <- predict(fit, type="response")
    sigma <- summary(fit)$sigma
    y.sim <- rnorm(n=length(y), mean=mu, sd=sigma)

  }

  return(y.sim)

}


#' Gaussian MX-knockoffs for continuous variables
#'
#' @param X data.frame (or tibble) with "numeric" columns only. The number of columns, ncol(X) needs to be > 2.
#'
#' @details \code{knockoffs_mx} performs MX knockoff simulation.
#'
#' @return Second-order multivariate Gaussian knockoff copy of X
#' @export
#'
#' @examples
#' #' library(knockofftools)
#'
#' set.seed(1)
#'
#' X <- generate_X(n=100, p=6, p_b=0, cov_type="cov_equi", rho=0.5)
#'
#' Xk <- knockoffs_mx(X)
knockoffs_mx <- function(X) {

  check_design(X, method="mx")

  knockoffs <- data.frame(knockoff::create.second_order(as.matrix(X)))

  return(knockoffs)

}

#' Knockoffs for general covariate data frames
#'
#' @param X data.frame (or tibble) with "numeric" and "factor" columns only. The number of columns, ncol(X) needs to be > 2.
#' @param method what type of knockoffs to calculate. Defaults to sequential knockoffs, method="seq", which works for
#' both numeric and factor variables. The other options are method="sparseseq", which is the sparse sequential knockoff generation algorithm, and method="mx", which only works if all columns of X are continuous.
#'
#' @return knockoff copy of X
#' @export
#'
#' @examples
#' library(knockofftools)
#'
#' set.seed(1)
#'
#' X <- generate_X(n=100, p=6, p_b=2, cov_type="cov_equi", rho=0.5)
#'
#' # sequential knockoff:
#' Xk <- knockoff(X)
#'
#' X <- generate_X(n=100, p=6, p_b=0, cov_type="cov_equi", rho=0.5)
#'
#' # MX-knockoff:
#' Xk <- knockoff(X, method="mx")
knockoff <- function(X, method="seq") {

  if (method=="seq") {
    Xk <- knockoffs_seq(X)
  }
  if (method=="mx") {
    Xk <- knockoffs_mx(X)
  }
  if (method=="sparseseq"){
    Xk <- knockoffs_sparse_seq(X)
  }

  return(Xk)

}
