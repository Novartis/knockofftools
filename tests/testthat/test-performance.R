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

test_that("fdp and tpp works", {
  expect_equal(eval_fdp(selected=c(1,2,3,5), negatives=5:10), 0.25)
  expect_equal(eval_tpp(selected=c(1,2,3,5), positives=1:4), 0.75)
})
