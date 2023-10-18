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

test_that("plot.variable.selections works as expected", {

  set.seed(1)

  # Simulate 8 Gaussian covariate predictors and 2 binary factors:
  X <- generate_X(n=100, p=10, p_b=2, cov_type="cov_equi", rho=0.2)

  # create linear predictor with first 5 beta-coefficients = 1 (all other zero)
  lp <- generate_lp(X, p_nn = 5, a=1)

  # Gaussian

  # Simulate response from a linear model y = lp + epsilon, where epsilon ~ N(0,1):
  y <- lp + rnorm(100)

  # Calculate M independent knockoff feature statistics:
  W <- knockoff.statistics(y=y, X=X, type="regression", M=10)

  S <-variable.selections(W, error.type = "pfer", level = 1)

  p <- suppressWarnings(plot(S))

  expect_equal(class(p), c("gg", "ggplot"))

  expect_error(plot.variable.selections(S$selected))

})
