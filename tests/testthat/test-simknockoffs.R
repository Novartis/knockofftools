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

test_that("simulation of knockoffs works", {

  set.seed(1)

  .gen.kock.mean <- function(setup) {
    with(setup, {
      X <- generate_X(n=n, p=p, p_b=p_b, cov_type=cov_type, rho=rho)

      Xk <- knockoff(X, method=method)

      rho <- mean(cor(Xk)[cor(Xk)!=1])

      expect_equal(rho, 0.5, tolerance=0.05)
    })
  }

  .gen.kock.mean(list(n=1000, p=6, p_b=0, cov_type="cov_equi", rho=0.5,
                      method="seq"))

  .gen.kock.mean(list(n=1000, p=6, p_b=0, cov_type="cov_equi", rho=0.5,
                      method="mx"))

  X <- generate_X(n=1000, p=6, p_b=2, cov_type="cov_equi", rho=0.5)

  expect_error(knockoff(X, method="mx"))

  Xk <- knockoff(X)

  expect_equal(lapply(X, is.factor), lapply(Xk, is.factor))
  expect_equal(cor(X[,unlist(lapply(X, is.numeric))]),
               cor(Xk[,unlist(lapply(Xk, is.numeric))]), tolerance=0.1)

})
