# Copyright 2023 Novartis Institutes for BioMedical Research Inc.
#
# Licensed under the MIT License (the "License");
# you may not use this file except in compliance with the License.
# You may obtain a copy of the License at
#
# https://www.mit.edu/~amini/LICENSE.md
#
# Unless required by applicable law or agreed to in writing, software
# distributed under the License is distributed on an "AS IS" BASIS,
# WITHOUT WARRANTIES OR CONDITIONS OF ANY KIND, either express or implied.
# See the License for the specific language governing permissions and
# limitations under the License.

test_that("Basic tests of knockoff.statistics and variable.selections", {

  set.seed(1)

  # Simulate 8 Gaussian covariate predictors and 2 (binary) factor predictors:
  X <- generate_X(n=100, p=10, p_b=2, cov_type="cov_equi", rho=0.2)

  # create linear predictor with first 5 beta-coefficients = 1 (all other zero)
  lp <- generate_lp(X, p_nn = 5, a=1)

  # Simulate response from a linear model y = lp + epsilon, where epsilon ~ N(0,1):
  y <- lp + rnorm(100)

  # Run the function with M=1
  W_M1 <- knockoff.statistics(y, X)

  expect_equal(class(W_M1), "data.frame")
  expect_equal(dim(W_M1), c(10, 1))
  expect_equal(names(W_M1), "W")

  .check.V <- function(W, M) {
    V <- variable.selections(W, error.type = "fdr", level = 0.2)

    expect_equal(class(V), c("variable.selections", "list"))
    expect_equal(class(V$selected), "data.frame")
    expect_equal(dim(V$selected), c(10, M))
    expect_equal(class(V$stable.variables), "character")

    V <- variable.selections(W, error.type = "fdr", level = 0)

    expect_equal(class(V$stable.variables), "character")
    expect_equal(length(V$stable.variables), 0)
  }


  .check.V(W_M1, 1)

  # Run the function with M=1
  W_M2 <- knockoff.statistics(y, X, M=2)

  .check.V(W_M2, 2)


  # Repeat the same checks for classification problems
  yb <- factor(rbinom(100, size=1, prob=exp(lp)/(1+exp(lp))))

  # Run the function with M=1
  W_M1 <- knockoff.statistics(yb, X, type = "classification")

  expect_equal(class(W_M1), "data.frame")
  expect_equal(dim(W_M1), c(10, 1))
  expect_equal(names(W_M1), "W")

  # Run the function with M=1
  W_M2 <- knockoff.statistics(yb, X, M=1, type = "classification")


  # Repeat the same checks for survival problems
  ys <- simulWeib(N=100, lambda0=0.01, rho=1, lp=lp)
  # Run the function with M=1
  W_M1 <- knockoff.statistics(ys, X, type = "survival")

  expect_equal(class(W_M1), "data.frame")
  expect_equal(dim(W_M1), c(10, 1))
  expect_equal(names(W_M1), "W")

  # Run the function with M=1
  W_M2 <- knockoff.statistics(ys, X, M=1, type = "survival")

})

test_that("Basic tests for stat_glmnet and stat_random_forest", {

  set.seed(1)

  # Simulate 10 Gaussian covariate predictors:
  X <- generate_X(n=100, p=10, p_b=2, cov_type="cov_equi", rho=0.2)

  #' # create linear predictor with first 5 beta-coefficients = 1 (all other zero)
  lp <- generate_lp(X, p_nn = 5, a=1)

  # Simulate response (g for gaussian, b for binary, and s for survival):
  yg <- lp + rnorm(100)
  yb <- factor(rbinom(100, size=1, prob=exp(lp)/(1+exp(lp))))
  ys <- simulWeib(N=100, lambda0=0.01, rho=1, lp=lp)

  # Calculate the knockoff copy of X (defaults to sequential knockoffs):
  X_k <- knockoff(X)

  # Calculate the knockoff statistic with stat_glmnet:
  Ws <- list(stat_glmnet(y=yg, X=X, X_k=X_k, type="regression"),
             stat_glmnet(y=yb, X=X, X_k=X_k, type="classification"),
             stat_glmnet(y=ys, X=X, X_k=X_k, type="survival"),
             stat_glmnet(y=yg, X=X, X_k=X_k,
                         X.fixed= data.frame(AGE=rnorm(nrow(X), mean=50, sd=10)),
                         type="regression"),
             stat_glmnet(y=yb, X=X, X_k=X_k,
                         X.fixed= data.frame(AGE=rnorm(nrow(X), mean=50, sd=10)),
                         type="classification"),
             stat_glmnet(y=ys, X=X, X_k=X_k,
                         X.fixed= data.frame(AGE=rnorm(nrow(X), mean=50, sd=10)),
                         type="survival"),
             stat_random_forest(X, X_k, yg, type = "regression"),
             stat_random_forest(X, X_k, yb, type = "classification"),
             stat_random_forest(X, X_k, ys, type = "survival"))

  .check.W <- function(W) {
    expect_equal(class(W), "data.frame")
    expect_equal(dim(W), c(10,1))
  }

  invisible(sapply(Ws, .check.W))

})


test_that("Basic tests for stat_predictive_glmnet and stat_predictive_causal_forest", {

  set.seed(1)

  # Simulate 10 Gaussian covariate predictors:
  X <- generate_X(n=100, p=10, p_b=2, cov_type="cov_equi", rho=0.2)

  #' # create linear predictor with first 5 beta-coefficients = 1 (all other zero)
  lp <- generate_lp(X, p_nn = 5, a=1)
  trt = sample(c(1,0), nrow(X), replace=TRUE)

  lp.pred = lp + 1*trt*( as.integer(X[,6]>0) + as.integer(X[,7]>0))

  # Simulate response (g for gaussian, b for binary, and s for survival):
  yg <- lp.pred + rnorm(100)
  yb <- factor(rbinom(100, size=1, prob=exp(lp.pred)/(1+exp(lp.pred))))
  ys <- simulWeib(N=100, lambda0=0.01, rho=1, lp=lp.pred)

  # Calculate the knockoff copy of X (defaults to sequential knockoffs):
  X_k <- knockoff(X)

  # Calculate the knockoff statistic with the tested filters
  Ws <- list(stat_predictive_glmnet(y=yg, X=X, X_k=X_k, trt=trt, type="regression"),
             stat_predictive_glmnet(y=yb, X=X, X_k=X_k, trt=trt, type="classification"),
             stat_predictive_glmnet(y=ys, X=X, X_k=X_k, trt=trt, type="survival"),
             stat_predictive_glmnet(y=yg, X=X, X_k=X_k, trt=trt,
                                    X.fixed= data.frame(AGE=rnorm(nrow(X), mean=50, sd=10)),
                                    type="regression"),
             stat_predictive_glmnet(y=yb, X=X, X_k=X_k, trt=trt,
                                    X.fixed= data.frame(AGE=rnorm(nrow(X), mean=50, sd=10)),
                                    type="classification"),
             stat_predictive_glmnet(y=ys, X=X, X_k=X_k, trt=trt,
                                    X.fixed= data.frame(AGE=rnorm(nrow(X), mean=50, sd=10)),
                                    type="survival"),
             stat_predictive_causal_forest(X, X_k, yg, trt=trt, type = "regression"),
             stat_predictive_causal_forest(X, X_k, yb, trt=trt, type = "classification"),
             stat_predictive_causal_forest(X, X_k, ys, trt=trt, type = "survival"))

  .check.W <- function(W) {
    expect_equal(class(W), "data.frame")
    expect_equal(dim(W), c(10,1))
  }

  invisible(sapply(Ws, .check.W))

})


test_that("Test different prognostic knockoff filters", {
  .check.WandV <- function(setup) {
    with(setup, {
      W <- knockoff.statistics(y=y, X=X, type=type, M=M, statistic=statistic,
                               trt=trt)

      expect_equal(class(W), "data.frame")
      expect_equal(dim(W), c(10, M))

      # Threshold the variable selections:
      V <- variable.selections(W, error.type = error.type, k=k, level = level)

      expect_equal(class(V), c("variable.selections", "list"))
    })
  }

  set.seed(0)

  # Simulate 10 Gaussian covariate predictors:
  X <- generate_X(n=100, p=10, p_b=2, cov_type="cov_equi", rho=0.2)

  #' # create linear predictor with first 5 beta-coefficients = 1 (all other zero)
  lp <- generate_lp(X, p_nn = 5, a=1)

  # Simulate response (g for gaussian, b for binary, and s for survival):
  y_g <- lp + rnorm(nrow(X))
  y_b <- factor(rbinom(nrow(X), size=1, prob=exp(lp)/(1+exp(lp))))
  y_surv <- simulWeib(N=nrow(X), lambda0=0.01, rho=1, lp=lp)

  expect_equal(class(y_surv), "Surv")

  # Test stat_glmnet for regression
  set.seed(1)
  .check.WandV(list(y=y_g, X=X, type="regression", M=5,
                    statistic="stat_glmnet", trt=NULL,
                    error.type = "pfer",k = NULL, level = 2))

  # Test stat_glmnet for classification
  set.seed(2)
  suppressWarnings({
    .check.WandV(list(y=y_b, X=X, type="classification", M=5,
                      statistic="stat_glmnet", trt=NULL,
                      error.type = "pfer",k = NULL, level = 2))
  })

  # Test stat_glmnet for survival
  set.seed(3)
  .check.WandV(list(y=y_surv, X=X, type="survival", M=5,
                    statistic="stat_glmnet", trt=NULL,
                    error.type = "kfwer", k=1,level=0.5))

  # Test stat_random_forest for regression
  set.seed(4)
  .check.WandV(list(y=y_g, X=X, type="regression", M=5,
                    statistic="stat_random_forest", trt=NULL,
                    error.type = "fdr",k = NULL, level = 0.2))

  # Test stat_random_forest for classification
  set.seed(5)
  .check.WandV(list(y=y_b, X=X, type="classification", M=5,
                    statistic="stat_random_forest", trt=NULL,
                    error.type = "fdr",k = NULL, level = 0.25))

  # Test stat_random_forest for survival
  set.seed(6)
  .check.WandV(list(y=y_surv, X=X, type="survival", M=5,
                    statistic="stat_random_forest", trt=NULL,
                    error.type = "pfer",level=1))
})


test_that("Test different predictive knockoff filters", {
  .check.WandV <- function(setup) {
    with(setup, {
      W <- knockoff.statistics(y=y, X=X, type=type, M=M, statistic=statistic,
                               trt=trt)

      expect_equal(class(W), "data.frame")
      expect_equal(dim(W), c(10, M))

      # Threshold the variable selections:
      V <- variable.selections(W, error.type = error.type, k=k, level = level)

      expect_equal(class(V), c("variable.selections", "list"))
    })
  }
  set.seed(0)

  # Simulate 10 Gaussian covariate predictors:
  X <- generate_X(n=100, p=10, p_b=2, cov_type="cov_equi", rho=0.2)

  # create prognostic part by a linear predictor with first 5 beta-coefficients = 1 (all other zero)
  lp_prog <- generate_lp(X, p_nn = 5, a=1)

  # Generate a binary treatment variable
  trt = sample(c(1,0), nrow(X), replace=TRUE)

  # create the final linear predictor introducing predictive part
  lp.pred = lp_prog + 1*trt*( as.integer(X[,6]>0) + as.integer(X[,7]>0))

  # Simulate response (g for gaussian, b for binary, and s for survival):
  y_g <- lp.pred + rnorm(nrow(X))
  y_b <- factor(rbinom(nrow(X), size=1, prob=exp(lp.pred)/(1+exp(lp.pred))))
  y_surv <- simulWeib(N=nrow(X), lambda0=0.01, rho=1, lp=lp.pred)

  # Test stat_predictive_glmnet for regression
  set.seed(1)
  .check.WandV(list(y=y_g, X=X, type="regression", M=5,
                    statistic="stat_predictive_glmnet", trt=trt,
                    error.type = "pfer", k=NULL,level=2))

  # Test stat_predictive_glmnet for classification
  set.seed(2)
  suppressWarnings({
    .check.WandV(list(y=y_b, X=X, type="classification", M=5,
                      statistic="stat_predictive_glmnet", trt=trt,
                      error.type = "pfer", k=NULL,level=2))
  })
  # Test stat_predictive_glmnet for survival
  set.seed(3)
  .check.WandV(list(y=y_surv, X=X, type="survival", M=5,
                    statistic="stat_predictive_glmnet", trt=trt,
                    error.type = "pfer", k=NULL,level=2))

  # Test stat_predictive_causal_forest for regression
  set.seed(4)
  .check.WandV(list(y=y_g, X=X, type="regression", M=5,
                    statistic="stat_predictive_causal_forest", trt=trt,
                    error.type = "pfer", k=NULL,level=2))

  # Test stat_predictive_causal_forest for classification
  set.seed(5)
  .check.WandV(list(y=y_b, X=X, type="classification", M=5,
                    statistic="stat_predictive_causal_forest", trt=trt,
                    error.type = "pfer", k=NULL,level=2))

  # Test stat_predictive_causal_forest for survival
  set.seed(6)
  suppressWarnings({
    .check.WandV(list(y=y_surv, X=X, type="survival", M=5,
                      statistic="stat_predictive_causal_forest", trt=trt,
                      error.type = "pfer", k=NULL, level=2))
  })
})


test_that("Test permutation based causal forest filters", {
  .check.WandV <- function(setup) {
    with(setup, {
      W <- knockoff.statistics(y=y, X=X, type=type, M=M, statistic=statistic,
                               trt=trt, permutations = permutations )

      expect_equal(class(W), "data.frame")
      expect_equal(dim(W), c(10, M))

      # Threshold the variable selections:
      V <- variable.selections(W, error.type = error.type, k=k, level = level)

      expect_equal(class(V), c("variable.selections", "list"))
    })
  }
  set.seed(0)

  # Simulate 10 Gaussian covariate predictors:
  X <- generate_X(n=250, p=10, p_b=2, cov_type="cov_equi", rho=0.2)

  # create prognostic part by a linear predictor with first 5 beta-coefficients = 1 (all other zero)
  lp_prog <- generate_lp(X, p_nn = 5, a=1)

  # Generate a binary treatment variable
  trt = sample(c(1,0), nrow(X), replace=TRUE)

  # create the final linear predictor introducing predictive part
  lp.pred = lp_prog + 1*trt*( as.integer(X[,6]>0) + as.integer(X[,7]==0))

  # Simulate response (g for gaussian and b for binary):
  y_g <- lp.pred + rnorm(nrow(X))
  y_b <- factor(rbinom(nrow(X), size=1, prob=exp(lp.pred)/(1+exp(lp.pred))))

  # Test stat_predictive_causal_forest for regression
  set.seed(1)
  .check.WandV(list(y=y_g, X=X, type="regression", M=1,
                    statistic="stat_predictive_causal_forest", trt=trt,
                    error.type = "pfer", k=NULL,
                    level=3,  permutations  = 5))

  # Test stat_predictive_causal_forest for classification
  set.seed(2)
  .check.WandV(list(y=y_b, X=X, type="classification", M=1,
                    statistic="stat_predictive_causal_forest", trt=trt,
                    error.type = "pfer", k=NULL,
                    level=2,    permutations  = 5))

})


test_that("Simple test of multi_select", {

  S <- matrix(c(1, 1, 0, 0, 0,
                1, 1, 0, 0, 0,
                0, 0, 0, 0, 0), ncol=3)

  expect_equal(multi_select(S), 1:2)

  S2 <- matrix(c(1, 1, 0, 0, 0,
                 1, 1, 1, 0, 0,
                 0, 0, 1, 0, 0), ncol=3)

  expect_equal(multi_select(S2), 3)

  S3 <- matrix(c(1, 1, 0, 0, 0,
                 1, 1, 0, 0, 0,
                 1, 1, 1, 0, 0,
                 1, 1, 1, 0, 0,
                 0, 0, 1, 0, 0), ncol=5)

  expect_equal(multi_select(S3), 1:2)

  S4 <- matrix(c(1, 1, 0, 0, 0,
                 1, 1, 0, 0, 0,
                 1, 1, 1, 0, 0,
                 1, 1, 1, 0, 0,
                 1, 1, 1, 0, 0), ncol=5)

  expect_equal(multi_select(S4), 1:3)

  S5 <- matrix(c(1, 1, 0, 0, 0,
                 1, 1, 0, 0, 0,
                 1, 1, 0, 0, 0,
                 1, 1, 1, 0, 0,
                 1, 1, 1, 0, 0,
                 1, 1, 1, 0, 0,
                 0, 1, 1, 0, 0), ncol=7)

  expect_equal(multi_select(S5), 1:2)

  S6 <- matrix(c(1, 1, 1, 1, 0,
                 1, 0, 1, 0, 0,
                 1, 1, 0, 0, 0,
                 1, 0, 1, 1, 0,
                 1, 1, 0, 1, 0,
                 1, 0, 1, 0, 1,
                 1, 1, 0, 0, 1), ncol=7)

  expect_equal(multi_select(S6), 1:2)

  S7 <- matrix(c(1, 1, 1, 1, 1,
                 1, 0, 1, 0, 0,
                 1, 1, 0, 0, 0,
                 1, 0, 1, 1, 0,
                 1, 1, 0, 1, 1,
                 1, 0, 1, 1, 1,
                 1, 1, 0, 0, 1), ncol=7)

  expect_equal(multi_select(S7), 1:2)

  S8 <- matrix(c(1, 1, 1, 1, 1,
                 1, 0, 1, 0, 1,
                 0, 0, 0, 0, 0,
                 1, 1, 1, 1, 0,
                 1, 1, 0, 1, 1,
                 1, 0, 1, 1, 0,
                 1, 1, 0, 0, 1), ncol=7)

  expect_equal(multi_select(S8), 1)

})


test_that("Expected errors and warnings", {

  set.seed(1)

  # Simulate 10 Gaussian covariate predictors:
  X <- generate_X(n=100, p=10, p_b=2, cov_type="cov_equi", rho=0.2)

  #' # create linear predictor with first 5 beta-coefficients = 1 (all other zero)
  lp <- generate_lp(X, p_nn = 5, a=1)

  # Simulate response (g for gaussian, b for binary, and s for survival):
  yg <- lp + rnorm(100)
  yb <- rbinom(100, size=1, prob=exp(lp)/(1+exp(lp)))
  ys <- simulWeib(N=100, lambda0=0.01, rho=1, lp=lp)

  expect_error(knockoff.statistics(yg, X, type="classification"))
  expect_error(knockoff.statistics(yg, X, type="survival"))
  expect_error(knockoff.statistics(factor(yb), X, type="regression"))
  expect_error(knockoff.statistics(factor(yb), X, type="survival"))
  expect_error(knockoff.statistics(yb, X, type="classification"))
  expect_error(knockoff.statistics(ys, X, type="regression"))
  expect_error(knockoff.statistics(ys, X, type="classification"))

  # Error, only accepts X as data.frame or tibble:
  expect_error(knockoff.statistics(yg, as.matrix(X)))

  # If try to supply binary numeric covariates then non-normality warning:
  X_with_numeric_factors <- dplyr::mutate_if(X, is.factor, as.numeric)
  expect_warning(knockoff.statistics(yg, X_with_numeric_factors),
                 "Some of the numeric columns of X have suspiciously few distinct values: n_distinct <= 30. Those columns should perhaps not be treated as continuous variables. Please review carefully and read the documentation about the gcm parameter of the knockoff.statistics function.")

  # Make one of the numerical variables skewed:
  X_with_7th_covariate_skewed <- X
  X_with_7th_covariate_skewed[,7] <- log(log(X_with_7th_covariate_skewed[,7]+4))
  expect_warning(knockoff.statistics(yg, X_with_7th_covariate_skewed, gcm = FALSE), "Some of the numeric input covariates may have ties and/or may not be normally distributed. This could affect the quality of corresponding knockoffs since they are sampled from a Gaussian distribution. X7 had normality rejected by Kolmogorov-Smirnov test. Please consider applying a normalizing transformation on these variables if needed.")

  # Create ties for first covariate:
  X_with_ties <- X
  X_with_ties[2,1] <- X_with_ties[1,1]
  expect_warning(knockoff.statistics(yg, X_with_ties, gcm = FALSE),
                 "Some of the numeric input covariates may have ties and/or may not be normally distributed. This could affect the quality of corresponding knockoffs since they are sampled from a Gaussian distribution. X1 had ties \\(but did not reject normality\\). Please consider applying a normalizing transformation on these variables if needed.")

  # What happens if there is only 2 covariates:
  set.seed(1)

  X <- generate_X(n=100, p=2, p_b=0, cov_type="cov_equi", rho=0.2)

  # Error when trying to compute knockoffs for a 2-column design matrix:
  expect_error(knockoff(X), "X should have ncol\\(X\\) > 2")

  # create linear predictor with first 5 beta-coefficients = 1 (all other zero)
  lp <- generate_lp(X, p_nn = 1, a=1)

  # Simulate response from a linear model y = lp + epsilon, where epsilon ~ N(0,1):
  y <- lp + rnorm(100)

  # Expect same error when calling knockoff.statistics:
  expect_error(knockoff.statistics(y, X), "X should have ncol\\(X\\) > 2")


})

# The following test tests the relevant gcm (normal score transformation)
# from the knockoff.statistics function:
test_that("Check if normal-score transformation is consistent", {

  set.seed(1)

  # Simulate 8 Gaussian covariate predictors and 2 (binary) factor predictors:
  X <- generate_X(n=100, p=10, p_b=0, cov_type="cov_equi", rho=0.2)

  mean.cor <- mean(cor(X)[cor(X)!=1])

  # Skew data heavily and then use the ns.transform code from knockoff.statistics function:
  X <- dplyr::mutate_all(X, function(x) exp(exp(x)))
  X <- dplyr::mutate_if(X, is.numeric, ns.transform)

  mean.cor.nstransformed <- mean(cor(X)[cor(X)!=1])

  expect_equal(round(mean.cor,2), round(mean.cor.nstransformed,2))

})
