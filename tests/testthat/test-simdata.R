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

test_that("generate_X and generate_y work", {

  set.seed(1)

  X <- generate_X(n=1000, p=10, p_b=0, cov_type="cov_diag", rho=0.5)

  expect_equal(class(X), "data.frame")
  expect_equal(dim(X), c(1000,10))
  expect_equal(sum(unlist(lapply(X, is.numeric))), 10)
  expect_equal(mean(cor(X)[cor(X)!=1]), 0, tolerance=0.05)

  X <- generate_X(n=1000, p=10, p_b=0, cov_type="cov_ar1", rho=0.5)

  expect_equal(class(X), "data.frame")
  expect_equal(dim(X), c(1000,10))
  expect_equal(sum(unlist(lapply(X, is.numeric))), 10)
  expect_equal(sum(abs(cor(X)[cor(X)!=1]-0.5) < 0.05) > 10, TRUE)
  expect_equal(sum(abs(cor(X)[cor(X)!=1]-0.5^2) < 0.05) > 10, TRUE)
  expect_equal(sum(abs(cor(X)[cor(X)!=1]-0.5^3) < 0.05) > 10, TRUE)

  X <- generate_X(n=1000, p=10, p_b=0, cov_type="cov_equi", rho=0.5)

  expect_equal(class(X), "data.frame")
  expect_equal(dim(X), c(1000,10))
  expect_equal(sum(unlist(lapply(X, is.numeric))), 10)
  expect_equal(mean(cor(X)[cor(X)!=1]), 0.5, tolerance=0.05)

  X <- generate_X(n=1000, p=10, p_b=10, cov_type="cov_equi", rho=0.5)

  expect_equal(class(X), "data.frame")
  expect_equal(dim(X), c(1000,10))
  expect_equal(sum(unlist(lapply(X, is.factor))), 10)

  X <- generate_X(n=1000, p=10, p_b=7, cov_type="cov_equi", rho=0.5)

  expect_equal(class(X), "data.frame")
  expect_equal(dim(X), c(1000,10))
  expect_equal(sum(unlist(lapply(X, is.numeric))), 3)
  expect_equal(sum(unlist(lapply(X, is.factor))), 7)

  y <- generate_y(X, p_nn=3, a=10)

  fm <- lm(y ~ ., data=cbind(y, X))

  expect_equal(as.numeric(summary(fm)$coefficients[2:4,4]), c(0, 0, 0), tolerance=0.0001)
  expect_equal(as.numeric(coefficients(fm)[2:4]), c(10, 10, 10), tolerance=0.1)
  expect_equal(as.numeric(coefficients(fm)[-c(2:4)]), rep(0, 8), tolerance=0.5)

})

test_that("simulWeib is sensible",{

  set.seed(1)

  X <- generate_X(n=10000, p=10, p_b=0, cov_type="cov_equi", rho=0.2)

  # create linear predictor with first 5 beta-coefficients = 1 (all other zero)
  lp <- generate_lp(X, p_nn = 5, a=10)

  # Simulate from Weibull hazard with with baseline hazard h0(t) = lambda*rho*t^(rho-1)
  # and linear predictor, whose first 3 coefficients are non-zero:
  y <- simulWeib(N=nrow(X), lambda0=0.01, rho=1, lp=lp)

  expect_equal(class(y), "Surv")

  # Fit cox regression:
  cx <- survival::coxph(y ~ ., data=cbind(y, X))

  expect_equal(as.numeric(coefficients(cx)[1:5]), rep(0.5, 5), tolerance=0.1)
  expect_equal(as.numeric(coefficients(cx)[6:10]), rep(0, 5), tolerance=0.1)

})
