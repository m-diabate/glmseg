## Submission summary

This is the first submission of **glmseg** to CRAN.

glmseg detects an unknown number of change-points in regression models (Gaussian
with homoscedastic or heteroscedastic noise, logistic, Poisson, and Weibull
accelerated failure time / survival), with covariates, via the max-EM algorithm.
The method is described in Diabate, Nuel and Bouaziz,
"Change-point detection in regression models via the max-EM algorithm",
Computational Statistics & Data Analysis (2025),
<doi:10.1016/j.csda.2025.108278>.

## R CMD check results

Local (macOS, R 4.2.1): 0 errors | 0 warnings | 1 note

* checking CRAN incoming feasibility ... NOTE
  Maintainer: 'Modibo Diabate <modibo.diabate.pro@gmail.com>'
  New submission
  Possibly misspelled words in DESCRIPTION:
    Bouaziz, Diabate, Nuel, initializations

This is the expected "New submission" note. The flagged words are the authors'
surnames (Bouaziz, Diabate, Nuel) and the correctly-spelled technical term
"initializations".

(A local-only note, "unable to verify current time", comes from the offline
build machine having no reliable network clock; it is not a package issue and
does not appear on CRAN's or win-builder's check machines.)

## Test environments

* local: macOS, R 4.2.1 -- 0 errors, 0 warnings
* win-builder, R-release and R-devel -- 0 errors, 0 warnings, 1 note
  (New submission)

## Reverse dependencies

There are currently no reverse dependencies (new package).

## Notes for the CRAN team

* The suggested package 'gfpop' is used only by the optional comparison helper
  gfpop_gauss(), guarded with requireNamespace(); the package works fully
  without it.
* Examples for exported functions that rely on suggested packages ('survival',
  'gfpop') are wrapped in if (requireNamespace(...)).
