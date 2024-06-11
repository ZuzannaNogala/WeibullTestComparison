
<!-- README.md is generated from README.Rmd. Please edit that file -->

# WeibullTestComparison

<!-- badges: start -->
<!-- badges: end -->

The repository contains R code used in my bachelor’s thesis, which
compares new tests proposed in work [New classes of tests for the
Weibull distribution using Stein’s method in the presence of random
right censoring
(2022)](https://link.springer.com/article/10.1007/s00180-021-01178-0) by
E. Bothma, J. S. Allison and I. J. H. Visagie with previously used test
like Kolmogorov-Smirnov, Cramér von Mises, Liao and Shimokawa and Krit.
However, in my work, I compared these tests on a broader range of
alternatives divided by the behavior of their hazard rate.

The repository contains the following files:

- `newton_alg.R`: Implements the Newton-Raphson algorithm to estimate
  the parameters of the Weibull distribution - $\lambda$ (scale
  parameter) and $\theta$ (shape parameter).
- `test_form.R`: Generates the values of the aforementioned tests
  without Krit test.
- `Quantiles_generator_uncens.R`: Contains all the necessary methods to
  generate quantiles for full samples.
- `powers_uncens.R`:Uses a `foreach` loop to generate the power of tests
  on a full sample, where the alternatives were the same as in the work
  mentioned above to check if the results are similar.
- `powers_uncens_HR.R` Uses a `foreach` loop to generate the power of
  tests on a full sample, with alternatives divided by the type of
  hazard rate.
- `powers_cens_HR.R`: Uses the bootstrap method to compute quantiles for
  censored data (from a uniform distribution) and a `foreach` loop to
  generate the power with alternatives divided by the type of hazard
  rate.

The goal is to understand the performance of these tests, particularly
in terms of power, when the important factor is the type of hazard rate.

This work was conducted in Polish in 2024 while I was a first-cycle
student at Wrocław University, Institute of Mathematics and Informatics.
