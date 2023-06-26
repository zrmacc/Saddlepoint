
# Saddlepoint Approximation

Zachary R. McCaw <br> Updated: 2023-06-26

## Description

This package performs the linear model score test with p-value
calculation via the empirical saddlepoint approximation.

## Installation

``` r
devtools::install_github(repo = "zrmacc/Saddlepoint")
```

## Usage

### Example data

Simulate example data for `n = 1e3` subjects where the proportion of
phenotypic variation explained by genotype is 0% and the proportion by
the covariate is 20%.

``` r
set.seed(100)
data <- Saddlepoint::DGP(n = 1e3, pve_g = 0.0, pve_x = 0.2)
```

### Fit the null model

Estimate the null model which excludes genotype. Extract the mean `mu`
and residual variance `sigma2` under the null.

``` r
fit <- stats::lm(y ~ x, data = data)
mu <- fit$fitted.values
sigma2 <- stats::sigma(fit)^2
```

### Perform score testing

- The required arguments are genotype `g`, phenotype `y`, and the null
  model mean `mu` and residual variance `sigma2`.
- The covariates `x` should either be supplied directly to the score
  test (as below), or else regressed from `g` prior to performing the
  score test.
- The output includes a flag `used_sp` of whether saddlepoint
  approximation was applied. If `FALSE`, the p-value was calculated via
  normal approximation. By default, normal approximation is applied when
  the magnitude of the score statistic is less than 1.96 (the 97.5th
  percentile of the normal distribution).

``` r
result <- Saddlepoint::ScoreTest(
  g = data$g,
  mu = mu,
  sigma2 = sigma2,
  x = data$x,
  y = data$y
)
show(result)
```

    ##   score_stat          p used_sp
    ## 1  -2.017887 0.04340702    TRUE
