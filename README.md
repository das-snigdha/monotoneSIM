Bayesian Estimation of Monotone Single Index Models
================

## Description and Intended Functionality

This package intends to perform a Markov Chain Monte Carlo algorithm to
estimate a monotone Single index models of the form:
*y* = *g*(*x*<sup>⊤</sup>*β*) + *ϵ*

where *y*: response, *x*: p-dimensional predictors, *β*: p-dimensional
coefficient vector, *ϵ*: independent N(0, *σ*<sub>*ϵ*</sub><sup>2</sup>)
random errors, *g*: unknown monotone link function. To construct a prior
on *g*, we consider the basis expansion:
*g*(*x*) = *Σ*<sub>*l* = 0</sub><sup>*L*</sup> *ξ*<sub>*l*</sub> *ψ*<sub>*l*</sub>(*x*),
where
*ψ*<sub>*l*</sub>(*x*) = ∫<sub> − 1</sub><sup>*x*</sup>*h*<sub>*l*</sub>(*t*) d*t*
, *h*<sub>*l*</sub> is a B-Spline basis function of order 2. Then, *g*
is monotonically increasing if and only if *ξ*<sub>*l*</sub> ≥ 0 ∀ *l*,
thereby enforcing the monotonicity constraint in an equivalent way.

The `monotoneSIM` function rescales the covariates and *β* so that
\|*x*<sup>⊤</sup>*β*\| ≤ 1. Then an MCMC algorithm is performed to
generate samples from the conditional posteriors of the unknown
parameters, *ξ*, *β*, *σ*<sub>*ϵ*</sub><sup>2</sup>. The obtained sample
of *β* is back-scaled to get *β* corresponding to the original
covariates. The function also returns the estimated monotone function
*g* for a grid of x values based on *ξ* obtained from the algorithm.

## Installation

To install this package from Github, run the following in your R
console:

``` r
devtools::install_github("das-snigdha/monotoneSIM")
```

## To Do

The following tasks shall be performed during the remaining of the
semester:

-   Perform Compatibility Checks on user supplied input.
-   Select equispaced `knots`,
    *u*<sub>0</sub>, *u*<sub>1</sub>, …, *u*<sub>*L*</sub> between -1
    and 1, if not supplied.
-   Select grid of x values, `grid.x` between -1 and 1 of length
    `size.grid.x`, if not supplied
-   Perform scaling of covariate matrix `X` and supplied `beta.init`.
-   Perform back-scaling on obtained `beta`, to get coefficients
    corresponding to unscaled X.
-   Create a vignette to demonstrate the usage, if time permits.
