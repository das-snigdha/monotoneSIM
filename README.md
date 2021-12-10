Bayesian Estimation of Monotone Single Index Models
================
Snigdha Das

## Overview

The package `monotoneSIM` performs Bayesian estimation using B-Spline
basis approximation of a Single Index Model where the unknown link
function is assumed to be monotonically increasing. Consider a Single
index models of the form:
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

Note: Here *g*(*x*) = 0, when *x* ≤  − 1 and *g*(*x*) = *g*(1), when
*x* ≥ 1. The matrix of predictors is scaled so that every row has
euclidean norm  ≤ 1 and euclidean norm of *β* is taken to be 1 so that
\|*x*<sup>⊤</sup>*β*\| ≤ 1.

## Functionality

The package provides two functions.

-   The main function `monotoneSIM` performs a Markov Chain Monte Carlo
    algorithm to generate samples from the conditional posterior
    distribution of unknown parameters in a monotone Single Index
    Model - the unknown parameters being B-Spline basis coefficients
    that approximate the unknown link function, the single index
    parameter and the error variance of the model.

-   The other function `monotoneFIT` calculates fitted values of the
    response using posterior mean of the single index parameter after
    generating an MCMC sample using the function `monotoneSIM`. It also
    returns the estimated link function of the model based on basis
    coefficients for a grid of supplied values.

## Installation

To install this package from Github, run the following in your R
console:

``` r
devtools::install_github("das-snigdha/monotoneSIM")
```

## An Example Demonstrating Usage

This package shall be implemented on the following simulated data.

We take two continuous variables and one dichotomous attribute as
predictors.

``` r
n = 100; p = 3;
X = matrix(rnorm(n*(p-1)), nrow = n, ncol = (p-1))
X = cbind(X, rbinom(n, 1, 0.5))
```

Let the true monotone link function be :
*g*(*x*) = 5\[*F*((*x* + 1)/2) − *F*(0)\],  − 1 ≤ *x* ≤ 1.

``` r
# True monotone link function
true.g = function(x){
  y = (x+1)/2
  5*(pnorm(y, mean=0.5, sd=0.1) - pnorm(0, mean = 0.5, sd = 0.1))
}
plot(true.g, xlim = c(-1,1))
```

Clearly g(.) is monotone increasing and *g*( − 1) = 0.

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
