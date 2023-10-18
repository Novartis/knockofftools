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

test_that("check_if_continuous works", {

  set.seed(1)
  X <- generate_X(n=1000, p=10, p_b=2, cov_type="cov_equi", rho=0.2)
  expect_equal(check_if_continuous(X), NULL)

  X$X1 <- sample(1:31, nrow(X), replace=TRUE)
  expect_equal(check_if_continuous(X), NULL)

  X$X1 <- sample(1:30, nrow(X), replace=TRUE)
  expect_warning(check_if_continuous(X))

})
