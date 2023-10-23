#' Internal check of whether class of response is compatible with family
#'
#' Do not call this function on its own
#'
#' @param y response variable
#' @param type "regression", "classification", or "survival"
#'
#' @return this function doesn't return anything
#' @export
#'
#' @keywords internal
check_family <- function(y, type) {

  if (!((type=="regression" & is.numeric(y) & class(y)!="Surv") | (type=="classification" & is.factor(y) & length(unique(y))==2) | (type=="survival" & class(y)=="Surv"))) {
    stop("One of the three must hold for the input: \n 1) type='regression' and class(y) = 'numeric', or \n 2) type='classification', class(y) = factor and length(unique(y)) = 2, or \n 3) type='survival' and class(y) = 'Surv' from library(survival)")
  }

}

#' Internal check of whether input data frame (or tibble) is of the right format
#'
#' Do not call this function on its own
#'
#' @param X data frame or tibble
#' @param method character string, either "seq" or "mx"
#'
#' @return this function doesn't return anything
#' @export
#'
#' @keywords internal
check_design <- function(X, method="seq", check.dim=TRUE) {

  if(!("data.frame" %in% class(X))) {
    stop(paste0(deparse(substitute(X)), " should be either a data.frame or tibble"))
  }

  if(check.dim & ncol(X)<=2) {
    stop(paste0(deparse(substitute(X)), " should have ncol(X) > 2"))
  }

  if(method=="seq" & sum(!unlist(lapply(X, class)) %in% c("factor", "numeric")) > 0) {
    stop(paste0(deparse(substitute(X)), " should only contain columns of class 'numeric' or 'factor'"))
  }

  if(method=="mx" & sum(!unlist(lapply(X, class)) %in% c("numeric")) > 0) {
    stop(paste0(deparse(substitute(X)), " should only contain columns of class 'numeric'"))
  }

}

#' Internal check of normality of the numeric input covariates
#'
#' Do not call this function on its own
#'
#' @param X data frame or tibble
#'
#' @return this function doesn't return anything
#' @export
#'
#' @keywords internal
check_normality <- function(X) {

  X_numeric <- dplyr::select_if(X, is.numeric)

  is.distinct <- unlist(lapply(X_numeric, dplyr::n_distinct))==nrow(X_numeric)

  p.values <- unlist(lapply(X_numeric, function(x) suppressWarnings(ks.test(x, y="pnorm", mean=mean(x), sd=sd(x))$p.value)))

  is.normal <- p.values >= 0.05

  if (!is.null(is.normal) & sum(!is.normal) > 0 | sum(!is.distinct) > 0) {

    warning.message <- paste0("Some of the numeric input covariates may have ties and/or may not be normally distributed. This could affect the quality of corresponding knockoffs since they are sampled from a Gaussian distribution. ")

    if (!is.null(is.normal) & sum(!is.normal) > 0) {
      warning.message <- paste0(warning.message, paste(names(which(!is.normal)), collapse=", "), " had normality rejected by Kolmogorov-Smirnov test. ")
    }

    if (sum(!is.distinct & is.normal) > 0) {
      warning.message <- paste0(warning.message, paste(names(which(!is.distinct & is.normal)), collapse=", "), " had ties (but did not reject normality). ")
    }

    warning.message <- paste0(warning.message, "Please consider applying a normalizing transformation on these variables if needed.")

    warning(warning.message)

  }

}

#' Select variables based on (heuristic) mode of multiple variable selections
#'
#' Do not call this function on its own
#'
#' @param S list of variable selection indices
#' @param p number of variables. Each element of the list of selection indices should be a subset of 1:p.
#' @param trim trimming probability threshold. A sensible default is \code{trim=0.5}.
#'
#' @return a single "most frequent" variable selection among the multiple selections in S.
#' @export
#'
#' @keywords internal
find_single_optimal_variable_set <- function(S, p, trim=0.5) {

  candidates <- 1:p

  countSpaces <- function(s) { sapply(gregexpr(" ", s), function(p) { sum(p>=0) } ) }

  Nknockoff <- length(S)

  var.freq.table <- table(factor(unlist(S), levels=candidates))

  which.remove <- as.numeric(which(var.freq.table < trim*Nknockoff))

  trimmed.selected <- lapply(S, function(x) paste(setdiff(x, which.remove), collapse=" "))

  model.freq.table <- table(unlist(trimmed.selected))

  best.trimmed.selected <- names(which(model.freq.table==max(model.freq.table)))

  # Resolve ties by choosing most parsimonious model:
  best.single.selected <- best.trimmed.selected[which.min(countSpaces(best.trimmed.selected))]

  selected <- as.integer(unlist(strsplit(best.single.selected, " ")))

  return(selected)

}


#' Internal function called within the stat_glmnet.
#'
#' Do not call this function on its own. Fits cross-validated glmnet model with fixed effect.
#'
#'
#' @param X.fixed a data.frame (or tibble) with "numeric" and "factor" columns corresponding to covariates or terms that should be treated as fixed effects in the model.
#' @param X original data.frame (or tibble) with "numeric" and "factor" columns only. The number of columns, ncol(X) needs to be > 2.
#' @param y response vector with \code{length(y) = nrow(X)}. Accepts "numeric" (family="gaussian") or binary "factor" (family="binomial"). Can also be a survival object of class Surv
#' as obtained from y = survival::Surv(time, status).
#' @param family should be "gaussian" if y is numeric, "binomial" if y is a binary factor variable or "cox" if y is a survival object.
#' @param nlambda length of lambda penalty sequence
#' @param penalty.factor Separate penalty factors, passed to the glmnet::cv.glmnet function, can be applied to each coefficient. This is a number that multiplies lambda to allow differential shrinkage. Can be 0 for some variables, which implies no shrinkage, and that variable is always included in the model.
#' @param ... other parameters passed to glmnet::cv.glmnet
#'
#' @return coefficients of both the original and knockoff variables
#' @export
#'
#' @keywords internal
cv_coeffs_glmnet_with_fixed_effect <- function(X_fixed, X, y, family, nlambda=500, penalty.factor, ...) {

  n = nrow(X)
  p = ncol(X)

  # Standardize variables
  if (!is.null(X_fixed)) {
    X <- scale(cbind(X, X_fixed))
  } else {
    X = scale(X)
  }

  if (!methods::hasArg(lambda) ) {
    if( identical(family, "gaussian") ) {
      if(!is.numeric(y)) {
        stop('Input y must be numeric.')
      }
      # Unless a lambda sequence is provided by the user, generate it
      lambda_max = max(abs(t(X) %*% y)) / n
      lambda_min = lambda_max / 2e3
      k = (0:(nlambda-1)) / nlambda
      lambda = lambda_max * (lambda_min/lambda_max)^k
    }
    else {
      lambda = NULL
    }
  }

  cv.glmnet.fit <- glmnet::cv.glmnet(x=X, y=y, lambda=lambda, standardize=F, standardize.response=F,
                                     family=family, penalty.factor = penalty.factor, ...)

  # If family == "cox" then there is no intercept term in the model:
  if (family == "cox") {
    coefs <- coef(cv.glmnet.fit, s = "lambda.min")[1:p]
  } else {
    coefs <- coef(cv.glmnet.fit, s = "lambda.min")[2:(p + 1)]
  }

  return(coefs)

}


#' Internal function called to return the importance scores from random forest
#'
#' @param X original data.frame with "numeric" and "factor" columns only.
#' @param y response vector with \code{length(y) = nrow(X)}. Accepts "numeric" (family="gaussian") or binary "factor" (family="binomial"). Can also be a survival object of class Surv
#' as obtained from y = survival::Surv(time, status).
#' @param type should be "regression" if y is numeric, "classification" if y is a binary factor variable or "survival" if y is a survival object.
#' @param ...
#'
#' @return importance scores
#' @export
#'
#' @keywords internal
random_forest_importance_scores <- function(X, y, trt, type = "regression", ...){
  # make the column names unique
  colnames(X) = make.unique(colnames(X))


  if (type=="survival") {
    df = data.frame(time=y[,1], status=y[,2], X)
    survival_formula <- as.formula(paste0("Surv(time, status) ~ ", paste(colnames(X), collapse="+")))
    rfFit = randomForestSRC::rfsrc(survival_formula, data=df, importance = TRUE, mtry = dim(X)[2])
    importance_scores = randomForestSRC::vimp.rfsrc(rfFit)$importance
  } else {
    df = data.frame(y, X)
    rfFit = randomForestSRC::rfsrc(y ~ ., data=df, importance=TRUE,  mtry = dim(X)[2])
    importance_scores = randomForestSRC::vimp.rfsrc(rfFit)$importance
  }


  return(importance_scores)
}


#' Internal function called to return the importance scores from causal forest
#'
#' @param X original data.frame with "numeric" and "factor" columns only.
#' @param y response vector with \code{length(y) = nrow(X)}. Accepts "numeric" (family="gaussian") or binary "factor" (family="binomial"). Can also be a survival object of class Surv
#' as obtained from y = survival::Surv(time, status).
#' @param trt binary treatment indicator
#' @param type should be "regression" if y is numeric, "classification" if y is a binary factor variable or "survival" if y is a survival object.
#' @param shuffle boolean variable, that takes the value FALSE if we do not want the target y to be shuffled, and TRUE otherwise
#' @param ...
#'
#' @return importance scores
#' @export
#'
#' @keywords internal
causal_forest_importance_scores <- function(X, y, trt, type = "regression", shuffle = FALSE,  ...){

  # make the column names unique
  colnames(X) = make.unique(colnames(X))
  X_transformed = model.matrix( ~ .+0, data=X, contrasts.arg = lapply(data.frame(X[,sapply(data.frame(X), is.factor)]),contrasts, contrasts = FALSE))

  if(type == "regression" | type ==  "classification"){
    if (shuffle == TRUE) {
      y = y[sample(dim(X)[1])]
    }
    c.forest = grf::causal_forest(X = X_transformed, Y = as.numeric(y), W = trt, mtry = dim(X_transformed)[2],  num.threads = 1, ...)
  }
  else if (type == "survival"){
    if (shuffle == TRUE) {
      y = y[sample(dim(X)[1]),]
    }
    funArgs <-  list(...)
    if (is.null(funArgs$horizon)){
      c.forest = grf::causal_survival_forest(X = X_transformed, Y = y[,1], W = trt, D = y[,2], mtry = dim(X_transformed)[2], horizon = min(max(y[trt==1,1]), max(y[trt==0,1])), num.threads = 1, ...)
    } else{
      c.forest = grf::causal_survival_forest(X = X_transformed, Y = y[,1], W = trt, D = y[,2], mtry = dim(X_transformed)[2], num.threads = 1,  ...)

    }
  }
  Z = grf::variable_importance(c.forest)

  # The columns of the model matrix X correspond to these original variables (stored in vars):
  indices_vars <- attributes(X_transformed)$assign
  # logical matrix that maps columns (of model.matrix) to variable indices:
  cols_to_vars <- outer(indices_vars, 1:max(indices_vars), function(x,y) x==y)

  # if j-th variable is factor with K levels then the variable importanze Z_j is:
  # sum(|beta_j,1|, ... , |beta_j,K-1|)

  importance_scores = apply(abs(c(Z)*cols_to_vars), 2, sum)

  return(importance_scores)
}




#' Ratio from Ren et al. (2021) used to derandomized knockoffs.
#'
#' Do not call this function on its own.
#'
#' @param M number of knockoff runs
#' @param eta the proportion of runs in which a feature
#' must be selected to be selected in the overall derandomized
#' procedure.
#'
#' @return returns pre-calculated ratios to be used in the multiple knockoff variable selection with error.type = "pfer" or "kfwer"
#' @export
#'
#' @keywords internal
find_ratio <- function(M, eta) {
  ratio = 1

  if (M > 5) {
    eta.ind <- which(round(eta, 2) == eta_list)
    M.ind <- which(M == M_list)
    ratio = M_eta_mat[M.ind, eta.ind]

  }

  return(ratio)
}



#' Normal score transformation function
#'
#' @param y a numeric vector representing the continuous variable
#'
#' @return a vector of length(y) with the normal-score transformed variable
#' @export
#'
#' @keywords internal
ns.transform <- function(y) {

  # Normal Q-Q plot:
  yt <- qqnorm(y, plot.it=FALSE)$x

  return(yt)

}

#' Heuristic check for whether a variable can be reasonably treated as continuous
#'
#' @param x a numeric variable vector
#'
#' @return a logical TRUE or FALSE depending on whether n_distinct(x) > 30
#' @export
#'
#' @keywords internal
check_if_continuous <- function(X) {
  `%>%` <- magrittr::`%>%`
  X_numeric <- dplyr::select_if(X, is.numeric)
  is.continuous <- sum(X_numeric %>% lapply(dplyr::n_distinct) %>% unlist() <= 30) > 0
  if (is.continuous) warning("Some of the numeric columns of X have suspiciously few distinct values: n_distinct <= 30. Those columns should perhaps not be treated as continuous variables. Please review carefully and read the documentation about the gcm parameter of the knockoff.statistics function.")
}


#' Estimate adjacency matrix using graphical LASSO (glasso)
#'
#' @param X data.frame (or tibble) of covariates
#'
#' @return adjacency matrix
#'
glasso_adjacency_matrix <- function(X){

  if (all(unlist(lapply(X, is.numeric)))){

    precision.matrix <- CVglasso::CVglasso(X, trace="none")$Omega # glasso::glasso(cov(X), rho=rho)$wi
    adjacency.matrix <- precision.matrix + t(precision.matrix)
    adjacency.matrix[adjacency.matrix !=0] = 1

  } else {

    # first, construct model matrix with dummy encoded factor variables
    X.dummy <- model.matrix(~., data=X)

    # The columns of the model matrix X correspond to these original variables (stored in vars):
    assignements <- attributes(X.dummy)$assign[-1]

    # Remove intercept of the model matrix
    X.dummy <- X.dummy[,-1]

    ## then construct its precision matrix
    precision.matrix.dummy <- CVglasso::CVglasso(X.dummy, trace="none")$Omega

    # This is done to guarantee symmetry in the adjacency matrix
    precision.matrix.dummy <- precision.matrix.dummy + t(precision.matrix.dummy)

    ## next construct adjacency matrix of original problem under the assumption that for each pair of nodes (A,B) involving at least one categorical variable, say A, there is an edge iff at least one of the dummy encoded features (of A) has an edge (to (any of the dummy encodings of ) B)
    adjacency.matrix.temp <- matrix(0, nrow=nrow(precision.matrix.dummy), ncol=ncol(X))
    ## combine entries by column
    for (i in c(1:ncol(X))){
      cols <- which(assignements==i)
      if (length(cols) == 1){
        adjacency.matrix.temp[,i] <- precision.matrix.dummy[,cols]
      } else {
        adjacency.matrix.temp[,i] = rowSums(precision.matrix.dummy[,cols])
      }
    }
    ## combine entries by row
    adjacency.matrix <- matrix(0, nrow=ncol(X), ncol=ncol(X))
    for (i in c(1:ncol(X))){
      rows <- which(assignements==i)
      if (length(rows) == 1){
        adjacency.matrix[i,] = adjacency.matrix.temp[rows,]
      } else {
        adjacency.matrix[i,] = colSums(adjacency.matrix.temp[rows,])
      }
    }

    ## convert to proper adjacency matrix (i.e. all entries either 1 or 0)
    adjacency.matrix[adjacency.matrix != 0] = 1

  }
  colnames(adjacency.matrix) <- names(X)
  return (adjacency.matrix)

}
