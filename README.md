# glmseg

[![License: GPL v3](https://img.shields.io/badge/License-GPLv3-blue.svg)](https://www.gnu.org/licenses/gpl-3.0)

**Change-point detection in regression models via the max-EM algorithm.**

`glmseg` detects an unknown number of change-points in *regression* models
(with covariates and non-Gaussian responses) using the max-EM algorithm, which
combines a constrained hidden Markov model with a classification-EM step.
Unlike exact dynamic-programming methods (PELT, FPOP, GFPOP), which operate on a
response vector, cannot include covariates, and (in the mean model) require
homoscedasticity, `glmseg` offers a single, consistent interface across a broad
class of regression models.

## Features

- Change-point detection for **Gaussian** (homoscedastic and heteroscedastic),
  **logistic**, **Poisson**, and **Weibull accelerated failure time**
  (right-censored survival) regression, all with **covariates**.
- Linear-time max-EM refinement, initialized by **Fused-Lasso** or **Binary
  Segmentation**.
- Selection of the number of segments by **BIC**, with an optional modified
  penalty (`"mbic"`) that curbs over-segmentation.
- An extensible **`family`** interface, full S3 methods, and comparison helpers
  (an exact brute-force search and a `gfpop` wrapper).

## Installation

From GitHub:

```r
# install.packages("remotes")
remotes::install_github("m-diabate/glmseg")
```

<!-- Once available on CRAN:
install.packages("glmseg")
-->

## Quick start

```r
library(glmseg)

set.seed(1)
y <- c(rnorm(200, 0, 1), rnorm(200, 8, 1), rnorm(200, -6, 1))

# Fixed number of segments (K = 3)
fit <- glmseg(y, K = 3, family = gaussian_homo(), init = "binseg")
fit$bp        #> 200 400
summary(fit)

# Choose the number of segments by BIC
select_k(y, K_max = 5, family = gaussian_homo(), penalty = "mbic")$K   #> 3
```

Other families follow the same interface, e.g.:

```r
# Logistic
yb <- c(rbinom(300, 1, plogis(-1.6)), rbinom(300, 1, plogis(1.6)))
glmseg(yb, K = 2, family = logistic_family())$bp

# Poisson
yc <- c(rpois(300, exp(0.4)), rpois(300, exp(2.2)))
glmseg(yc, K = 2, family = poisson_family())$bp

# Survival (Weibull AFT, right-censored), requires the 'survival' package
# survseg(time, status, K = 2)
```

See the vignette (`vignette("glmseg")`) for a tour of all families.

## Method

The max-EM algorithm is described in:

> Diabaté, M., Nuel, G., and Bouaziz, O. (2025).
> *Change-point detection in regression models via the max-EM algorithm.*
> Computational Statistics & Data Analysis.
> <https://doi.org/10.1016/j.csda.2025.108278>

## License

GPL-3. See [LICENSE](https://www.gnu.org/licenses/gpl-3.0).
