#' Knockoff (feature) statistics:
#'
#' This function calculates M >= 1 independent knockoff (feature) statistics (W) given input response vector y and covariate data.frame X.
#' The function first calculates M independent knockoff copies (Xk1, ..., XkM) of the covariate matrix (X) and then calculates the knockoff
#' feature statistics W1, ..., WM. By default each feature statistic is calculated via the parameter
#' statistic=knockofftools::stat_glmnet, but user may write and supply their own feature statistics functions
#' (e.g. random forest variable importance difference). The user may additionally supply fixed effects (X.fixed)
#' that should always be included in the underlying model (e.g. covariates to adjust for).
#'
#' If multiple knockoffs are desired (M > 1) then the method utilizes the clustermq package for parallel distribution of jobs on HPC scheduler. See \href{https://mschubert.github.io/clustermq/articles/userguide.html}{clustermq-userguide}
#' for further details on how to configure (differently from defaults) clustermq scheduler and batch templates.
#'
#' @param y response vector with \code{length(y) = nrow(X)}. Accepts "numeric", binary "factor", or survival ("Surv") object.
#' @param X data.frame (or tibble) with "numeric" and "factor" columns only. The number of columns, ncol(X) needs to be > 2.
#' @param type should be "regression" if y is numeric, "classification" if y is a binary factor variable or "survival" if y is a survival object.
#' @param M the number of independent knockoff feature statistics that should be calculated.
#' @param knockoff.method what type of knockoffs to calculate. Defaults to sequential knockoffs, knockoff.method="seq", with other options: knockoff.method="sparseseq" and knockoff.method="mx".
#' The "mx" method only works if all columns of the X matrix are continuous.
#' @param statistic knockoff feature statistic function, defaults to glmnet coefficient difference (statistic="stat_glmnet"; see ?stat_glmnet). Other options include statistic="stat_random_forest" (see ?stat_random_forest), statistic="stat_predictive_glmnet" (see ?stat_predictive_glmnet) or statistic="stat_predictive_causal_forest" (see ?stat_predictive_causal_forest).
#' @param trt binary treatment (factor) variable required if statistic involves a predictive knockoff filter (i.e. if statistic="stat_predictive_glmnet" or statistic="stat_predictive_causal_forest")
#' @param gcm logical indicator for whether a Gaussian Copula Model should be applied. Defaults to TRUE since the underlying knockoff generation mechanism for numeric variables is based on multivariate Gaussian variables.
#' When gcm=TRUE each numeric variable is normal score transformed resulting in marginal standard normal variables. The knockoff filter then acts in this transformed variable space. User is advised not to change this parameter unless he/she understands the consequences.
#' @param ... additional parameters passed to the "statistic" function (note that the knockoffs parameter X_k should not be entered by user; it is already calculated inside the knockoff.statistics function).
#'
#' @return data.frame with knockoff statistics W as column. The number of rows matches the number of columns (variables) of the data.frame X and the variable names are recorded in rownames(W).
#' @export
#'
#' @examples
#' library(knockofftools)
#'
#' set.seed(1)
#'
#' # Simulate 10 Gaussian covariate predictors:
#' X <- generate_X(n=100, p=10, p_b=0, cov_type="cov_equi", rho=0.2)
#'
#' # create linear predictor with first 5 beta-coefficients = 1 (all other zero)
#' lp <- generate_lp(X, p_nn = 5, a=1)
#'
#' # Gaussian
#'
#' # Simulate response from a linear model y = lp + epsilon, where epsilon ~ N(0,1):
#' y <- lp + rnorm(100)
#'
#' # Calculate M independent knockoff feature statistics:
#' W <- knockoff.statistics(y=y, X=X, type="regression", M=5)
#' print(variable.selections(W, error.type = "pfer", level = 2)$stable.variables)
#' \dontrun{
#' W <- knockoff.statistics(y=y, X=X, type="regression", M=5, statistic = "stat_random_forest")
#' print(variable.selections(W, error.type = "pfer", level = 2)$stable.variables)
#'
#' # Cox
#'
#' # Simulate from Weibull hazard with with baseline hazard h0(t) = lambda*rho*t^(rho-1)
#' # and linear predictor lp:
#' y <- simulWeib(N=nrow(X), lambda0=0.01, rho=1, lp=lp)
#'
#' W <- knockoff.statistics(y=y, X=X, type="survival", M=5)
#' print(variable.selections(W, error.type = "pfer", level = 2)$stable.variables)
#'
#' W <- knockoff.statistics(y=y, X=X, type="survival", M=5, statistic = "stat_random_forest")
#' print(variable.selections(W, error.type = "pfer", level = 2)$stable.variables)
#'
#' # Check quickly predictive filters
#' # Generate a binary treatment variable
#' trt = sample(c(1,0), nrow(X), replace=TRUE)
#' lp.pred = lp + 1*trt*( as.integer(X[,6]>0) + as.integer(X[,7]>0))
#'
#' # Gaussian
#'
#' y <- lp.pred + rnorm(nrow(X))
#'
#' W <- knockoff.statistics(y=y, X=X, type="regression", statistic = "stat_predictive_glmnet", trt=trt, M=5)
#' print(variable.selections(W, error.type = "pfer", level = 2)$stable.variables)
#'
#' W <- knockoff.statistics(y=y, X=X, type="regression", statistic = "stat_predictive_causal_forest", trt=trt, M=5)
#' print(variable.selections(W, error.type = "pfer", level = 2)$stable.variables)
#'
#' # Cox
#'
#' # Simulate from Weibull hazard with with baseline hazard h0(t) = lambda*rho*t^(rho-1)
#' # and linear predictor lp:
#' y <- simulWeib(N=nrow(X), lambda0=0.01, rho=1, lp=lp.pred)
#'
#' W <- knockoff.statistics(y=y, X=X, type="survival", statistic = "stat_predictive_glmnet", trt=trt, M=5)
#' print(variable.selections(W, error.type = "pfer", level = 2)$stable.variables)
#'
#' W <- knockoff.statistics(y=y, X=X, type="survival",
#'                          statistic = "stat_predictive_causal_forest",
#'                          trt=trt, M=5)
#' print(variable.selections(W, error.type = "pfer", level = 2)$stable.variables)}
knockoff.statistics <- function(y, X, type="regression",
                                M = 1, knockoff.method = "seq",
                                statistic = "stat_glmnet", trt=NULL, gcm=TRUE,
                                ...) {

  check_design(X)

  # Check if the numeric columns of X can be regarded as "continuous"
  # (simple heuristic check whether the number of distinct values are above 30)
  check_if_continuous(X)

  # By default apply Gaussian Copula Model via normal score transformation on the numeric columns of X
  if (gcm) {
    X <- dplyr::mutate_if(X, is.numeric, ns.transform)
  } else{
    # If gcm is not used, check normality
    check_normality(X)
  }


  check_family(y, type)

  if (grepl("predictive", statistic) & is.null(trt)) {
    stop("statistic=", statistic,
         " is a predictive filter and requires the binary factor treatment",
         " variable 'trt' to be specified.")
  }

  if (!grepl("predictive", statistic) & !is.null(trt)) {
    warning("statistic=", statistic,
            " is not a predictive filter, hence the 'trt' variable",
            " will be ignored.")
  }

  if (M==1) {
      W <- .knockoff.statistics.single(y, X, type=type,
                                        knockoff.method = knockoff.method,
                                        statistic = statistic, trt=trt,
                                        ...)
  } else if (round(M)==M & M > 1) {
    W <- clustermq::Q(.knockoff.statistics.single,
                      type=rep(type, M),
                      const = list(y=y, X=X, knockoff.method=knockoff.method, statistic=statistic, trt=trt, ...),
                      pkgs = "knockofftools",
                      n_jobs=M,
                      log_worker = FALSE,
                      timeout = 300)
    W <- as.data.frame(W)
    names(W) <- paste0("W", 1:ncol(W))
  } else {
    stop("M needs to be an integer >= 1")
  }

  return(W)

}

#' Knockoff (feature) statistics for a single knockoff
#'
#' Do not call this function on its own. Use \code{knockoff.statistics} instead.
#'
#' @param y response vector with \code{length(y) = nrow(X)}. Accepts "numeric", binary "factor", or survival ("Surv") object.
#' @param X data.frame (or tibble) with "numeric" and "factor" columns only. The number of columns, ncol(X) needs to be > 2.
#' @param type should be "regression" if y is numeric, "classification" if y is a binary factor variable or "survival" if y is a survival object.
#' @param M the number of independent knockoff feature statistics that should be calculated.
#' @param knockoff.method what type of knockoffs to calculate. Defaults to sequential knockoffs, knockoff.method="seq", but other options are "sparseseq" and "mx". The "mx" option only works if all columns of X are continuous.
#' @param statistic knockoff feature statistic function, defaults to glmnet coefficient difference (statistic="stat_glmnet"; see ?stat_glmnet). Other options include statistic="stat_random_forest" (see ?stat_random_forest), statistic="stat_predictive_glmnet" (see ?stat_predictive_glmnet) or statistic="stat_predictive_causal_forest" (see ?stat_predictive_causal_forest).
#' @param trt binary treatment (factor) variable required if statistic involves a predictive knockoff filter (i.e. if statistic="stat_predictive_glmnet" or statistic="stat_predictive_causal_forest")
#' @param ... additional parameters passed to the "statistic" function (note that the knockoffs parameter X_k should not be entered by user; it is already calculated inside the knockoff.statistics function).
#'
#' @return data.frame with a single knockoff statistics W as column.
#' @export
#'
#' @keywords internal
.knockoff.statistics.single <- function(y, X, type="regression",
                                        knockoff.method = "seq",
                                        statistic = "stat_glmnet", trt=NULL,
                                        ...) {

  # Calculate the knockoff copy of X (defaults to sequential knockoffs):
  X_k <- knockoff(X, method=knockoff.method)
  # Calculate the knockoff statistic from a single knockoff copy:
  statistic <- get(statistic)
  W <- statistic(y=y, X=X, X_k=X_k, type=type, trt=trt, ...)

  return(W)

}

#' Knockoff (feature) statistics: Absolute elastic-net coefficient differences between original and knockoff variables
#'
#' This function follows mostly the implementation of knockoff::glmnet.stat_coefdiff. The input data.frames (X, X_k) and X.fixed (if supplied)
#' are first converted to design matrices (with the function model.matrix). This means that if
#' the input features contain factor variables then associated dummy variable are determined by the model.matrix contrasts
#' (defaults to indicator dummy variables with a reference level). There is then a call to glmnet::cv.glmnet where the input is
#' y and x = cbind(X, X_k, X.fixed) and penalty is only applied to cbind(X, X_k). If user wishes to also penalize the parameters
#' of X.fixed then an additional penalty.fixed parameter can be adjusted accordingly.
#'
#' If there are factor covariates with multiple levels among columns of X then there will be more columns in model.matrix
#' than in the corresponding data.frame (both for original X and its knockoff X_k). In this case, let W_j be the difference
#' between the two maximum absolute signals of coefficients of model.matrix associated with covariate j. I.e. if j-th variable
#' is factor with K levels then W_j is: max(|beta_j1|, ... , |beta_j,K-1|) - max(|beta.tilde_j1|, ..., |beta.tilde_j,K-1|) where
#' (beta_j1, ..., beta_j,K-1) and (beta.tilde_j1, ..., beta.tilde_j,K-1) are the coefficients associated with dummy variables
#' of original j-th factor and its knockoff, respectively.
#'
#' @param X original data.frame (or tibble) with "numeric" and "factor" columns only. The number of columns, ncol(X) needs to be > 2.
#' @param X_k knockoff data.frame (or tibble) with "numeric" and "factor" columns only obtained e.g. by X_k = knockoff(X). The dimensions and column classes must match
#' those of the original X.
#' @param y response vector with \code{length(y) = nrow(X)}. Accepts "numeric" (type="regression") or binary "factor" (type="classification"). Can also be a survival object of class "Surv" (type="survival")
#' as obtained from y = survival::Surv(time, status).
#' @param type should be "regression" if y is numeric, "classification" if y is a binary factor variable or "survival" if y is a survival object.
#' @param X.fixed a data.frame (or tibble) with "numeric" and "factor" columns corresponding to covariates or terms that should be treated as fixed effects in the model.
#' @param penalty.fixed a numeric vector of length equal to number of columns of X.fixed indicating which fixed effects should be estimated with glmnet penalty and which not
#' (1 corresponds to covariates that should be penalized and 0 corresponds to covariates that are not penalized; if X.fixed is supplied, all elements of penalty.fixed are set to zero as default)
#' @param ... additional parameters passed to glmnet::cv.glmnet
#'
#' @return data.frame with knockoff statistics W as column. The number of rows matches the number of columns (variables) of the data.frame X and the variable names are recorded in rownames(W).
#' @export
#'
#' @examples
#' library(knockofftools)
#'
#' set.seed(1)
#'
#' # Simulate 10 Gaussian covariate predictors and 1 factor with 4 levels:
#' X <- generate_X(n=100, p=10, p_b=0, cov_type="cov_equi", rho=0.2)
#'
#' # Simulate response from a linear model y = X%*%beta + epsilon, where epsilon ~ N(0,1) with
#' # first 3 beta-coefficients = 1 (all other zero):
#' y <- (X$X1 + X$X2 + X$X3) + rnorm(100)
#'
#' # Calculate M independent knockoff feature statistics:
#' W <- knockoff.statistics(y=y, X=X)
stat_glmnet <- function(y, X, X_k, type = "regression", X.fixed=NULL, penalty.fixed = rep(0, length(X.fixed)), ...) {

  if (type=="regression") family <- "gaussian"
  if (type=="classification") family <- "binomial"
  if (type=="survival") family <- "cox"

  check_design(X); check_design(X_k); if(!is.null(X.fixed)) check_design(X.fixed, check.dim=FALSE)

  # Randomly swap columns of X and Xk
  swap = as.logical(rbinom(ncol(X), 1, 0.5))
  X.swap <- X; X.swap[,swap] <- X_k[,swap]
  Xk.swap <- X_k; Xk.swap[,swap] <- X[,swap]

  # Calculate model matrices of the data.frames X and X_k and X.fixed [if supplied] (recycle X, X_k, and X.fixed):
  X.swap <- model.matrix(~., data=X.swap)
  Xk.swap <- model.matrix(~., data=Xk.swap)

  # The columns of the model matrix X correspond to these original variables (stored in vars):
  indices_vars <- attributes(X.swap)$assign[-1]

  if (!is.null(X.fixed)) {
    X.fixed <- model.matrix(~., data=X.fixed)
    # The columns of the model matrix X.fixed correspond to these original variables:
    indices_vars_fixed <- attributes(X.fixed)$assign[-1]
    # Expand penalty.fixed to all columns of model.matrix (e.g. all levels will inherit the penalty of a given factor variable):
    penalty.fixed <- penalty.fixed[indices_vars_fixed]
  }

  # Remove intercepts of the model matrices (this -1 is importantly also in definition of indices_vars and indices_vars_fixed):
  X.swap <- X.swap[,-1]; Xk.swap <- Xk.swap[,-1]; X.fixed <- X.fixed[,-1]

  penalty.factor <- c(rep(1, ncol(X.swap) + ncol(Xk.swap)), penalty.fixed)

  # Compute statistics
  Z = cv_coeffs_glmnet_with_fixed_effect(X.fixed, cbind(X.swap, Xk.swap), y, family=family, penalty.factor = penalty.factor, ...)

  p = ncol(X.swap)
  orig = 1:p

  # If the columns of model.matrix map one-to-one to columns of data.frame then use regular solution from knockoff::stat.glmnet_coefdiff:
  if (p == ncol(X)) {
    W = abs(Z[orig]) - abs(Z[orig+p])
  } else {

    # if j-th variable is factor with K levels then W_j is calculated as follows:
    # W.vec_j <- c(|beta_j1|-|beta.tilde_j1|, ... , |beta_j,K-1|-|beta.tilde_j,K-1|)
    # W_j = W.vec_j[which.max(abs(W.vec_j))]

    # logical matrix that maps columns (of model.matrix) to variable indices:
    cols_to_vars <- outer(indices_vars, 1:max(indices_vars), function(x,y) x==y)

    W = apply(abs(Z[orig]*cols_to_vars) - abs(Z[orig+p]*cols_to_vars), 2, function(x) x[which.max(abs(x))])

  }

  # Correct for swapping of columns of X and Xk
  W = W * (1-2*swap)

  # Return a named vector (with variable names):
  W.dataframe <- data.frame(W = W)
  rownames(W.dataframe) <- names(X)

  return(W.dataframe)

}


#' Knockoff (feature) statistics: Random forest
#'
#' This function uses the ranger package to estimate permutation importance scores from random forest.
#'
#' If there are factor covariates with multiple levels among columns of X then there will be more columns in model.matrix
#' than in the corresponding data.frame (both for original X and its knockoff X_k). In this case, let W_j be the difference
#' between the two sums derived by the variable importance (VI) scores associated with covariate j. I.e. if j-th variable
#' is factor with K levels then W_j is: sum(|VI_j,1|, ... , |VI_j,K|) - sum(|VI_j1|, ..., |VI_j,K|).
#'
#' @param X original data.frame (or tibble) with "numeric" and "factor" columns only. The number of columns, ncol(X) needs to be > 2.
#' @param X_k knockoff data.frame (or tibble) with "numeric" and "factor" columns only obtained e.g. by X_k = knockoff(X). The dimensions and column classes must match
#' those of the original X.
#' @param y response vector with \code{length(y) = nrow(X)}. Accepts "numeric" (family="gaussian") or binary "factor" (family="binomial"). Can also be a survival object of class "Surv" (type="survival")
#' as obtained from y = survival::Surv(time, status).
#' @param type should be "regression" if y is numeric, "classification" if y is a binary factor variable or "survival" if y is a survival object.
#' @param ...
#'
#' @return data.frame with knockoff statistics W as column. The number of rows matches the number of columns (variables) of the data.frame X and the variable names are recorded in rownames(W).
#' @export
#'
#' @examples
#' library(knockofftools)
#'
#' set.seed(1)
#'
#' # Simulate 10 Gaussian covariate predictors and 1 factor with 4 levels:
#' X <- generate_X(n=500, p=10, p_b=0, cov_type="cov_diag", rho=0.2)
#' X$X11 <- factor(sample(c("A","B","C","D"), nrow(X), replace=TRUE))
#'
#' # Calculate the knockoff copy of X:
#' X_k <- knockoff(X)
#'
#' # create linear predictor with first 3 beta-coefficients = 1 (all other zero) and a treatment effect of size 1
#' lp <- (X$X1 + X$X2 + X$X3)
#'
#' # Gaussian
#'
#' # Simulate response from a linear model y = lp + epsilon, where epsilon ~ N(0,1):
#' y <- lp + rnorm(nrow(X))
#'
#' W <- stat_random_forest(X, X_k, y, type = "regression")
#'
#' # Cox
#'
#' # Simulate from Weibull hazard with with baseline hazard h0(t) = lambda*rho*t^(rho-1) and linear predictor lp:
#' y <- simulWeib(N=nrow(X), lambda0=0.01, rho=1, lp=lp)
#'
#' # Calculate  knockoff feature statistics:
#' W <- stat_random_forest(X, X_k, y, type = "survival")
#'
stat_random_forest <- function(X, X_k, y, type = "regression",  ...) {

  check_design(X); check_design(X_k)

  # Randomly swap columns of X and Xk
  swap = as.logical(rbinom(ncol(X), 1, 0.5))
  X.swap <- X; X.swap[,swap] <- X_k[,swap]
  Xk.swap <- X_k; Xk.swap[,swap] <- X[,swap]

  p = ncol(X.swap)
  orig = 1:p
  # Run in the original data
  var_split_imps= random_forest_importance_scores(cbind(X.swap,Xk.swap), y, type = type, ...)

  W = var_split_imps[orig] - var_split_imps[orig+p]


  # Correct for swapping of columns of X and Xk
  W = W * (1-2*swap)

  # Return a named vector (with variable names):
  W.dataframe <- data.frame(W = W)
  rownames(W.dataframe) <- names(X)

  return(W.dataframe)

}


#' Knockoff (feature) statistics that captues the predictive strength: Absolute coefficient differences between treatment original variables interaction terms and treatment knockoff variables interaction terms
#'
#' This function follows the implementation of the stat_glmnet function, but modifies it to focus on the coefficients of the interaction terms with the treatment.
#' A first version of this filter was suggested in Sechidis et al. (2021), but with this function we also allow the user to adjust to variables, or to
#' fix some variables as prognostic in the model.
#'
#' If there are factor covariates with multiple levels among columns of X then there will be more columns in model.matrix
#' than in the corresponding data.frame (both for original X and its knockoff X_k). In this case, let W_j be the difference
#' between the two maximum absolute signals of coefficients of model.matrix associated with covariate j. I.e. if j-th variable
#' is factor with K levels then W_j is: max(|beta_j1|, ... , |beta_j,K-1|) - max(|beta.tilde_j1|, ..., |beta.tilde_j,K-1|) where
#' (beta_j1, ..., beta_j,K-1) and (beta.tilde_j1, ..., beta.tilde_j,K-1) are the coefficients associated with dummy variables
#' of original j-th factor and its knockoff, respectively.
#'
#' @param X original data.frame (or tibble) with "numeric" and "factor" columns only. The number of columns, ncol(X) needs to be > 2.
#' @param X_k knockoff data.frame (or tibble) with "numeric" and "factor" columns only obtained e.g. by X_k = knockoff(X). The dimensions and column classes must match
#' those of the original X.
#' @param y response vector with \code{length(y) = nrow(X)}. Accepts "numeric" (type="regression") or binary "factor" (type="classification"). Can also be a survival object of class "Surv" (type="survival")
#' as obtained from y = survival::Surv(time, status).
#' @param trt a binary treatment indicator variable (should be numeric with 0/1 entries)
#' @param type should be "regression" if y is numeric, "classification" if y is a binary factor variable or "survival" if y is a survival object.
#' @param X.fixed a data.frame (or tibble) with "numeric" and "factor" columns corresponding to covariates or terms that should be treated as fixed effects in the model.
#' @param penalty.fixed a numeric vector of length equal to number of columns of X.fixed indicating which fixed effects should be estimated with glmnet penalty and which not
#' (1 corresponds to covariates that should be penalized and 0 corresponds to covariates that are not penalized; if X.fixed is supplied, all elements of penalty.fixed are set to zero as default)
#' @param fixed.prognostic a parameters that describes weather the user fix some prognostic variable in the model (TRUE) or not (FALSE)
#' @param ... additional parameters passed to glmnet::cv.glmnet
#'
#' @return data.frame with knockoff statistics W as column that capture the predictive strength of the variables. The number of rows matches the number of columns (variables) of the data.frame X and the variable names are recorded in rownames(W).
#' @export
#'
#' @examples
#' \dontrun{
#' library(knockofftools)
#'
#' set.seed(1)
#'
#' # Simulate 10 Gaussian covariate predictors and 1 factor with 4 levels:
#' X <- generate_X(n=500, p=10, p_b=0, cov_type="cov_diag", rho=0.2)
#' X$X11 <- factor(sample(c("A","B","C","D"), nrow(X), replace=TRUE))
#'
#'  # Calculate the knockoff copy of X:
#' X_k <- knockoff(X)
#'
#' # Generate a binary treatment variable
#' trt = sample(c(1,0), nrow(X), replace=TRUE)
#'
#' # Simulate a fixed "treatment" effect:
#' X.fixed <- data.frame(SEX = factor(sample(c("male", "female"), nrow(X), replace=TRUE)), trt = trt)
#' penalty.fixed = rep(0, length(X.fixed))
#'
#' # create linear predictor with first 3 beta-coefficients = 1 (all other zero) and a treatment effect of size 1
#' lp <- X.fixed$trt+ as.numeric(X.fixed$SEX == 'male') + (X$X1 + X$X2 + X$X3) + (X$X4 + 2*as.integer(X$X11=='A'))*trt
#'
#' # Gaussian
#'
#' # Simulate response from a linear model y = lp + epsilon, where epsilon ~ N(0,1):
#' y <- lp + rnorm(nrow(X))
#'
#' # We present three different ways of using that filter
#' # (a) the filter suggested in Sechidis et al. (2021) - where there is no fixed prognostic part
#' W <- stat_predictive_glmnet(y=y, trt=trt, X=X, X_k=X_k, type="regression", X.fixed=X.fixed, penalty.fixed = penalty.fixed, fixed.prognostic = FALSE)
#' # (b) a variation where the user can fix some prognostic variables in the model and allows the model to penalise them
#' X.fixed.prognostic <- data.frame( X1 = X$X1, X2 = X$X2, X3 = X$X3)
#' penalty.fixed.prognostic = rep(1, length(X.fixed.prognostic))
#' W <- stat_predictive_glmnet(y=y, trt=trt, X=X, X_k=X_k, type="regression", X.fixed=cbind(X.fixed, X.fixed.prognostic),  penalty.fixed = c(penalty.fixed,penalty.fixed.prognostic), fixed.prognostic = TRUE)
#'
#' # Cox
#'
#' # Simulate from Weibull hazard with with baseline hazard h0(t) = lambda*rho*t^(rho-1) and linear predictor lp:
#' y <- simulWeib(N=nrow(X), lambda0=0.01, rho=1, lp=lp)
#'
#' # Explore the three scenarios above with survival outcomes
#' # (a) the filter suggested in Sechidis et al. (2021) - where there is no fixed prognostic part
#' W <- stat_predictive_glmnet(y=y, trt=trt, X=X, X_k=X_k, type="survival", X.fixed=X.fixed, penalty.fixed = penalty.fixed, fixed.prognostic = FALSE)
#' # (b) a variation where the user can fix some prognostic variables in the model and allows the model to penalise them
#' X.fixed.prognostic <- data.frame( X1 = X$X1, X2 = X$X2, X3 = X$X3)
#' penalty.fixed.prognostic = rep(1, length(X.fixed.prognostic))
#' W <- stat_predictive_glmnet(y=y, trt=trt, X=X, X_k=X_k, type="survival", X.fixed=cbind(X.fixed, X.fixed.prognostic),  penalty.fixed = c(penalty.fixed,penalty.fixed.prognostic), fixed.prognostic = TRUE)
#'}
#' @details Sechidis, K., Kormaksson, M., & Ohlssen, D. (2021). Using knockoffs for controlled predictive biomarker identification. Statistics in Medicine, 40(25), 5453-5473.
stat_predictive_glmnet <- function(X, X_k, y, trt, type = "regression", X.fixed=NULL, penalty.fixed = rep(0, length(X.fixed)), fixed.prognostic = FALSE,...) {

  if (type=="regression") family <- "gaussian"
  if (type=="classification") family <- "binomial"
  if (type=="survival") family <- "cox"

  check_design(X); check_design(X_k); if(!is.null(X.fixed)) check_design(X.fixed, check.dim=FALSE)


  # Randomly swap columns of X and Xk
  swap = as.logical(rbinom(ncol(X), 1, 0.5))
  X.swap <- X; X.swap[,swap] <- X_k[,swap]
  Xk.swap <- X_k; Xk.swap[,swap] <- X[,swap]

  # Calculate model matrices of the data.frames X and X_k and X.fixed [if supplied] (recycle X, X_k, and X.fixed):
  X.swap <- model.matrix(~., data=X.swap)
  Xk.swap <- model.matrix(~., data=Xk.swap)

  # Combine the original Xs and their interactions T:X
  X.swap.interactions = matrix(0, dim(X.swap)[1], dim(X.swap)[2])
  for (col in 1:dim(X.swap)[2]){
    X.swap.interactions[,col]=  X.swap[,col]*trt
  }

  # Generate interactions with knockoffs the knockoffs Xs and their interactions T:X
  Xk.swap.interactions = matrix(0, dim(Xk.swap)[1], dim(Xk.swap)[2])
  for (col in 1:dim(Xk.swap)[2]){
    Xk.swap.interactions[,col] =  Xk.swap[,col]*trt
  }

  # The columns of the model matrix X correspond to these original variables (stored in vars):
  indices_vars <- attributes(X.swap)$assign[-1]

  if (!is.null(X.fixed)) {
    X.fixed <- model.matrix(~., data=X.fixed)
    # The columns of the model matrix X.fixed correspond to these original variables:
    indices_vars_fixed <- attributes(X.fixed)$assign[-1]
    # Expand penalty.fixed to all columns of model.matrix (e.g. all levels will inherit the penalty of a given factor variable):
    penalty.fixed <- penalty.fixed[indices_vars_fixed]
  }

  # Remove intercepts of the model matrices (this -1 is importantly also in definition of indices_vars and indices_vars_fixed):
  X.swap <- X.swap[,-1]; Xk.swap <- Xk.swap[,-1]; X.swap.interactions <- X.swap.interactions[,-1]; Xk.swap.interactions <- Xk.swap.interactions[,-1]; X.fixed <- X.fixed[,-1]

  # Compute statistics
  if (fixed.prognostic){
    penalty.factor <- c(rep(1, ncol(X.swap.interactions) + ncol(Xk.swap.interactions)), penalty.fixed)
    Z_predictive = cv_coeffs_glmnet_with_fixed_effect(X.fixed, cbind(X.swap.interactions, Xk.swap.interactions), y, family=family, penalty.factor = penalty.factor, ...)

    p = ncol(X.swap)
    orig = 1:p

    # If the columns of model.matrix map one-to-one to columns of data.frame then use regular solution from knockoff::stat.glmnet_coefdiff:
    if (p == ncol(X)) {
      W_predictive = abs(Z_predictive[orig]) - abs(Z_predictive[orig+p])

    } else {

      # if j-th variable is factor with K levels then W_j is:
      # max(|beta_j1|, ... , |beta_j,K-1|) - max(|beta.tilde_j1|, ..., |beta.tilde_j,K-1|)

      # logical matrix that maps columns (of model.matrix) to variable indices:
      cols_to_vars <- outer(indices_vars, 1:max(indices_vars), function(x,y) x==y)

      W_predictive = apply(abs(Z_predictive[orig]*cols_to_vars), 2, max) - apply(abs(Z_predictive[orig+p]*cols_to_vars), 2, max)

    }
  } else{
    penalty.factor <- c(rep(1, ncol(X.swap) + ncol(Xk.swap)+ ncol(X.swap.interactions) + ncol(Xk.swap.interactions)), penalty.fixed)
    Z = cv_coeffs_glmnet_with_fixed_effect(X.fixed, cbind(X.swap, Xk.swap, X.swap.interactions, Xk.swap.interactions), y, family=family, penalty.factor = penalty.factor, ...)

    p = ncol(X.swap)
    orig = 1:p

    Z_predictive = Z[(2*p+1):(4*p)]
    # If the columns of model.matrix map one-to-one to columns of data.frame then use regular solution from knockoff::stat.glmnet_coefdiff:
    if (p == ncol(X)) {
      W_predictive = abs(Z_predictive[orig]) - abs(Z_predictive[orig+p])

    } else {

      # if j-th variable is factor with K levels then W_j is:
      # max(|beta_j1|, ... , |beta_j,K-1|) - max(|beta.tilde_j1|, ..., |beta.tilde_j,K-1|)

      # logical matrix that maps columns (of model.matrix) to variable indices:
      cols_to_vars <- outer(indices_vars, 1:max(indices_vars), function(x,y) x==y)

      W_predictive = apply(abs(Z_predictive[orig]*cols_to_vars), 2, max) - apply(abs(Z_predictive[orig+p]*cols_to_vars), 2, max)

    }
  }
  # Correct for swapping of columns of X and Xk
  W_predictive = W_predictive * (1-2*swap)

  # Return a named vector (with variable names):
  W_predictive.dataframe <- data.frame(W = W_predictive)
  rownames(W_predictive.dataframe) <- names(X)


  return(W_predictive.dataframe)

}



#' Causal forest based knockoff (feature) statistics that captues the predictive strength: Difference from importance scores derived by causal forest
#'
#' This filter presented in Sechidis et al. (2021).
#'
#' If there are factor covariates with multiple levels among columns of X then there will be more columns in model.matrix
#' than in the corresponding data.frame (both for original X and its knockoff X_k). In this case, let W_j be the difference
#' between the two sums derived by the variable importance (VI) scores associated with covariate j. I.e. if j-th variable
#' is factor with K levels then W_j is: sum(|VI_j,1|, ... , |VI_j,K|) - sum(|VI_j1|, ..., |VI_j,K|).
#'
#' @param X original data.frame (or tibble) with "numeric" and "factor" columns only. The number of columns, ncol(X) needs to be > 2.
#' @param X_k knockoff data.frame (or tibble) with "numeric" and "factor" columns only obtained e.g. by X_k = knockoff(X). The dimensions and column classes must match
#' those of the original X.
#' @param y response vector with \code{length(y) = nrow(X)}. Accepts "numeric" (type="regression") or binary "factor" (type="classification"). Can also be a survival object of class "Surv" (type="survival")
#' as obtained from y = survival::Surv(time, status).
#' @param trt a binary treatment indicator variable (should be numeric with 0/1 entries)
#' @param type should be "regression" if y is numeric, "classification" if y is a binary factor variable or "survival" if y is a survival object.
#' @param permutations when it is not null, an integer that defines the number of permutations for the method suggested in O'Neil and Weeks (2018).
#' @param ... additional parameters passed to grf::causal_forest (for type = "regression" and "classification) and causal_survival_forest (for type = "survival")
#'
#' @return data.frame with knockoff statistics W as column that capture the predictive strength of the variables. The number of rows matches the number of columns (variables) of the data.frame X and the variable names are recorded in rownames(W).
#' @export
#'
#' @examples
#' library(knockofftools)
#'
#' set.seed(1)
#'
#' # Simulate 10 Gaussian covariate predictors and 1 factor with 4 levels:
#' X <- generate_X(n=500, p=10, p_b=0, cov_type="cov_diag", rho=0.2)
#' X$X11 <- factor(sample(c("A","B","C","D"), nrow(X), replace=TRUE))
#'
#'  # Calculate the knockoff copy of X:
#' X_k <- knockoff(X)
#'
#' # Generate a binary treatment variable
#' trt = sample(c(1,0), nrow(X), replace=TRUE)
#'
#' # Simulate a fixed "treatment" effect:
#' X.fixed <- data.frame(SEX = factor(sample(c("male", "female"), nrow(X), replace=TRUE)), trt = trt)
#' penalty.fixed = rep(0, length(X.fixed))

#' # create linear predictor with first 3 beta-coefficients = 1 (all other zero) and a treatment effect of size 1
#' lp <- X.fixed$trt+ as.numeric(X.fixed$SEX == 'male') + (X$X1 + X$X2 + X$X3) + (X$X4 + as.integer(X$X11=='A'))*trt
#'
#' # Gaussian
#'
#' # Simulate response from a linear model y = lp + epsilon, where epsilon ~ N(0,1):
#' y <- lp + rnorm(nrow(X))
#'
#' # Without permutations
#' W <- stat_predictive_causal_forest(X=X, X_k=X_k, y=y, trt=trt, type="regression")
#'
#' # With permutations
#' W <- stat_predictive_causal_forest( X=X, X_k=X_k, y=y, trt=trt, type="regression", permutations = 100)
#'
#' # Cox
#'
#' # Simulate from Weibull hazard with with baseline hazard h0(t) = lambda*rho*t^(rho-1)
#' # and linear predictor lp:
#' y <- simulWeib(N=nrow(X), lambda0=0.01, rho=1, lp=lp)
#'
#' # Without permutations
#' W <- stat_predictive_causal_forest(X=X, X_k=X_k, y=y, trt=trt, type="survival")
#'
#' # With permutations
#' W <- stat_predictive_causal_forest(X=X, X_k=X_k, y=y, trt=trt, type="survival", permutations = 100)
#'
#' @details Sechidis, K., Kormaksson, M., & Ohlssen, D. (2021). Using knockoffs for controlled predictive biomarker identification. Statistics in Medicine, 40(25), 5453-5473.
#' @details O'Neill, E., & Weeks, M. (2018). Causal tree estimation of heterogeneous household response to time-of-use electricity pricing schemes. arXiv preprint arXiv:1810.09179.
stat_predictive_causal_forest <- function(X, X_k, y, trt, type = "regression", permutations = NULL, ...) {
  check_design(X); check_design(X_k)


  # Randomly swap columns of X and Xk
  swap = as.logical(rbinom(ncol(X), 1, 0.5))
  X.swap <- X; X.swap[,swap] <- X_k[,swap]
  Xk.swap <- X_k; Xk.swap[,swap] <- X[,swap]

  p = ncol(X.swap)
  orig = 1:p
  # Run in the original data
  var_split_imps= causal_forest_importance_scores(cbind(X.swap,Xk.swap), y, trt, type = type,  shuffle = FALSE,  ...)

  if(is.null(permutations)){
    W = var_split_imps[orig] - var_split_imps[orig+p]
  } else{

    n_jobs <- 2
    perm_var_split_imps <- clustermq::Q(causal_forest_importance_scores,
                                        type = rep(type,permutations),
                                        const = list(X=cbind(X,X_k),
                                                     y=as.matrix(y),
                                                     trt = trt,
                                                     shuffle = TRUE, ...),
                                        pkgs = c("knockofftools","grf"),
                                        n_jobs=n_jobs)

    mat_perm_var_split_imps <- do.call(cbind, perm_var_split_imps)
    test_bind <- cbind(var_split_imps,mat_perm_var_split_imps)
    p_values <-  apply(test_bind, 1,
                       function(x) sum(x[1] < x[2:ncol(test_bind)])/permutations)

    W = -p_values[orig] + p_values[orig+p]
  }

  # Correct for swapping of columns of X and Xk
  W = W * (1-2*swap)

  # Return a named vector (with variable names):
  W_predictive.dataframe <- data.frame(W = W)
  rownames(W_predictive.dataframe) <- names(X)


  return(W_predictive.dataframe)

}





#' Controls the false discovery rate (FDR) given knockoff W-statistics.
#'
#' This method was introduced by Candes et al. (2018)
#'
#' @param W a vector of knockoff W-statistics (feature statistics).
#' @param level The nominal level at which to control the FDR.
#'
#' @return the selected variables
#' @export
#'
#' @examples
#' W <- rnorm(100)
#' selections_control_FDR(W, level=0.5)
#'
#' @details E. Candes, Y. Fan, L. Janson, & J. Lv, (2018). Panning for gold:‘model‐X’knockoffs for high dimensional controlled variable selection. Journal of the Royal Statistical Society: Series B (Statistical Methodology), 80(3), 551-577.
selections_control_FDR <- function(W, level) {
  W <- as.numeric(W)

  ko.thres <- knockoff::knockoff.threshold(W, fdr=level, offset=1)
  selected <- which(W >= ko.thres)

  return(selected)

}

#' Controls the per-familywise error rate (PFER) given knockoff W-statistics.
#'
#' This method was introduced by Janson and Su (2016), while we used the implementation of https://github.com/zhimeir/derandomized_knockoffs_paper/blob/master/R/pfer_filter.R
#'
#' @param W a vector of knockoff W-statistics (feature statistics)
#' @param level The nominal level at which to control the PFER.
#'
#' @return the selected variables
#' @export
#'
#' @examples
#' W <- rnorm(100)
#' selections_control_PFER(W, level=2)
#'
#' @details L. Janson, & W. Su, (2016). Familywise error rate control via knockoffs. Electronic Journal of Statistics, 10(1), 960-975.
selections_control_PFER <- function(W, level){
  v0 <- level
  W <- as.numeric(W)
  v <- floor(v0)+rbinom(1,1,v0-floor(v0))
  order_w <- order(abs(W),decreasing = TRUE)
  sorted_w <- W[order_w]
  negid <- which(sorted_w<0)
  S <- c()
  if(v>0){
    if(v> length(negid)){
      #if the total number of negitives is less than v, select all positives
      S <-  which(W>0)
    }else{
      TT <- negid[v]
      S <- which(sorted_w[1:TT]>0)
      S <- order_w[S]
    }
  } else {
    S <- integer(0)
  }
  return(S)

}




#' Controls the k-familywise error rate (k-FWER) given a vector of knockoff W-statistics.
#'
#' This method was introduced by Janson and Su (2016),, while we used the implementation https://github.com/zhimeir/derandomized_knockoffs_paper/blob/master/R/vanilla_fwer_filter.R
#'
#' @param W a vector of knockoff W-statistics (feature statistics)
#' @param level The nominal level at which to control the k-FWER.
#' @param k  a positive integer corresponding to k-FWER (multiple testing when one seeks to control at least k false discoveries)
#'
#' @return the selected variables
#' @export
#'
#' @examples
#' W <- rnorm(100)
#' selections_control_kFWER(W, level=0.1, k=5)
#'
#' @details L. Janson, & W. Su, (2016). Familywise error rate control via knockoffs. Electronic Journal of Statistics, 10(1), 960-975.
selections_control_kFWER <- function(W, level, k) {
  # if k = 0, we can't reject anything
  if (k == 0) {return(c())}


  v_list <- seq(1,20,by=1)
  p_list <- pnbinom(k,v_list,.5,lower.tail = FALSE)+dnbinom(k,v_list,.5)
  if(sum(p_list<=level)>0){
    v0 <- max(which(p_list<=level))
    prob <- (level-p_list[v0])/(p_list[v0+1] - p_list[v0])
  }else{
    v0 <- 0
    prob <- (level)/(p_list[1])
  }

  ## Run v-knockoffs
  v <- floor(v0)+rbinom(1,1,prob)

  S <- selections_control_PFER(W=W, level=v)
  return(S)
}





#' Knockoff variable selection: Select the variables by controlling a user-specified error rate
#'
#' This is the main function that performs the knockoff based variable selection using as input the knockoff statistics W.
#' In case of multiple knockoffs, ncol(W) > 1, the function performs variable selection for each knockoff and
#' additionally stabilizes the selections by combining their outcomes.
#'
#' @param W a data.frame of knockoff W-statistics (feature statistics); columns correspond to different knockoffs and rows correspond to the underlying variables. row.names(W) records the variable names.
#' @param level the nominal level that the user wants to control
#' @param error.type the error rate to control, at the moment "fdr", "pfer" and "kfwer"
#' @param k  a positive integer corresponding to k-FWER (multiple testing when one seeks to control at least k false discoveries), to be used only with error.type = 'kfwer'
#' @param thres threshold parameter for stabilizing the selections (eta parameter for derandomized knockoffs, trims parameter for multi_select). A natural choice is thres = 0.5.
#'
#' @details Knockoffs is a randomized procedure which relies on the construction of synthetic (knockoff) variables.
#' This function performs variable selection for multiple knockoffs and then stabilizes the selections by combining their outcomes.
#' When the pfer or kfwer error is controlled the derandomizing knockoffs is used, which was introduced by Ret et al. (2021) and provably controls this errors.
#' When the fdr is controlled the heuristic multiple selection algorithm is used, which was introduced by Kormaksson et al. (2021).
#'
#' @return an object of class "variable.selections" that is essentially a list with two elements: 1) $selections = (p x M) binary data.frame where rows correspond to variables, and cols correspond to different knockoffs; a value of 1 means
#' the given variable was selected for that particular knockoff simulation, 0 otherwise; 2) $stable.selection =  a character vector with the selected variables from stability selection (as described in Details). The second field is
#' only meaningful if user specifies multiple knockoffs (say M > 5). If M = 1 then the stable.selection simply returns the indicies of $selections that are equal to 1.
#' @export
#'
#' @seealso \code{\link[knockofftools]{plot.variable.selections}} for plotting an organized heatmap of the selections.
#'
#' @examples
#' library(knockofftools)
#'
#' set.seed(1)
#'
#' # Simulate 10 Gaussian covariate predictors:
#' X <- generate_X(n=100, p=10, p_b=0, cov_type="cov_equi", rho=0.2)
#'
#' # create linear predictor with first 5 beta-coefficients = 1 (all other zero)
#' lp <- generate_lp(X, p_nn = 5, a=1)
#'
#' # Gaussian
#'
#' # Simulate response from a linear model y = lp + epsilon, where epsilon ~ N(0,1):
#' y <- lp + rnorm(100)
#'
#' # Calculate M independent knockoff feature statistics:
#' W <- knockoff.statistics(y=y, X=X, type="regression", M=5)
#'
#' S = variable.selections(W, error.type = "pfer", level = 1)
#'
#' # selections under alternative error control:
#' S = variable.selections(W, error.type = "kfwer", k=1, level = 0.50)
#' S = variable.selections(W, error.type = "fdr", level = 0.5)
#'
#' @details Z. Ren, Y. Wei, & E. Candès, (2021). Derandomizing knockoffs. Journal of the American Statistical Association, 1-11.
#' @details M. Kormaksson, L. J. Kelly, X. Zhu, S. Haemmerle, L. Pricop, & D. Ohlssen (2021). Sequential knockoffs for continuous and categorical predictors: With application to a large psoriatic arthritis clinical trial pool. Statistics in Medicine, 40(14), 3313-3328.
variable.selections <- function(W, level = 0.20, error.type = "fdr", k = NULL, thres=0.50) {
  ## Check the type of error criterion
  if(error.type %in% c("fdr","pfer","kfwer") == 0) stop("The error criterion is not supported!")

  # Preprocessing
  p <- nrow(W)
  M <- ncol(W)
  error.type <- tolower(error.type)

  # Choose appropriate variable selection function (which.select) and in the case of "pfer" and "kfwer" adjust nomal level w.r.t. M and thres
  if (error.type == "pfer") {
    which.select <- selections_control_PFER
    ratio <-  find_ratio(M, thres)
    level = level/ratio
  }
  if (error.type == "kfwer") {
    which.select <- ifelse(M==1, selections_control_kFWER, selections_control_PFER)
    ratio <- find_ratio(M, thres)
    level <- ifelse(M==1, level, ifelse(k==1, level/ratio, 2*k*level/ratio))
  }
  if (error.type == "fdr") {
    which.select = selections_control_FDR
  }

  # Loop through W-statistics to generate the binary matrix of selections S
  S = matrix(0, p, M)
  for (i in 1:M) {
    if ((error.type == "kfwer")&(M==1)) {
      S[which.select(W[,i], level = level, k = k),i] <- 1
    } else{
      S[which.select(W[,i], level = level),i] <- 1
    }
  }

  # Perform the final selection
  if (error.type == 'pfer'| error.type == 'kfwer') {
    selected_variables = which(rowMeans(S)>thres)
  } else if (error.type == 'fdr'){
    selected_variables = multi_select(S = S, trim = thres)
  }

  S <- as.data.frame(S)
  rownames(S) <- rownames(W)
  names(S) <- gsub("W", "S", names(W))

  object <- list(selected = S, stable.variables = rownames(S)[selected_variables])

  class(object) <- c("variable.selections", class(object))

  return(object)

}



#' Select variables based on the heuristic multiple selection algorithm from Kormaksson et al. 'Sequential
#' knockoffs for continuous and categorical predictors: With application to a large psoriatic arthritis clinical
#' trial pool.' Statistics in Medicine. 2021;1–16.
#'
#' @param S the binary matrix of selections
#' @param trim trimming probability threshold. A sensible default is \code{trim=0.5}.
#'
#' @return a single "most frequent" variable selection among the multiple selections in S.
#' @export
#'
#' @examples
#' library(knockofftools)
#'
#' set.seed(1)
#'
#' p = 31
#' Nknockoff = 100
#' S <- matrix(sample(0:1,p*Nknockoff, replace=TRUE), p, Nknockoff)
#'
#' multi_select(S)
#' @details M. Kormaksson, L. J. Kelly, X. Zhu, S. Haemmerle, L. Pricop, & D. Ohlssen (2021). Sequential knockoffs for continuous and categorical predictors: With application to a large psoriatic arthritis clinical trial pool. Statistics in Medicine, 40(14), 3313-3328.
multi_select <- function(S, trim=0.5) {

  p <- nrow(S)
  Nknockoff <- ncol(S)
  candidates <- 1:p

  var.freq.table <- rowMeans(S)

  trim.seq <- unique(var.freq.table)
  trim.seq <- trim.seq[trim.seq >= trim]

  selection_list <- lapply(as.data.frame(S), function(x) which(x==1))

  if (length(trim.seq) > 0) {
    selected <- lapply(trim.seq, function(x) find_single_optimal_variable_set(selection_list, p, trim=x))
    selected <- selected[[which.max(unlist(lapply(selected, length)))]]
  } else {
    selected <- integer(0)
  }

  return(selected)

}



