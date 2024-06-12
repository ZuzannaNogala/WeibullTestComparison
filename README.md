
<!-- README.md is generated from README.Rmd. Please edit that file -->

# WeibullTestComparison

<!-- badges: start -->
<!-- badges: end -->

\[ENG\]

The repository contains R code used in my bachelor’s thesis, which
compares new tests proposed in work [New classes of tests for the
Weibull distribution using Stein’s method in the presence of random
right censoring
(2022)](https://link.springer.com/article/10.1007/s00180-021-01178-0) by
E. Bothma, J. S. Allison and I. J. H. Visagie with previously used test
like Kolmogorov-Smirnov, Cramér von Mises, [Liao and
Shimokawa](https://www.tandfonline.com/doi/abs/10.1080/00949659908811965)
and [Krit](http://www.numdam.org/item/JSFS_2014__155_3_135_0/). However,
in my work, I compared these tests on a broader range of alternatives
divided by the behavior of their hazard rate.

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

\[POL\]

Repozytorium zawiera kod napisany w R użyty w mojej pracy licencjackiej,
gdzie zostały porównane nowe testy zaproponowane w artykule [New classes
of tests for the Weibull distribution using Stein’s method in the
presence of random right censoring
(2022)](https://link.springer.com/article/10.1007/s00180-021-01178-0)
autorstwa E. Bothma, J. S. Allison i I. J. H. Visagie z wcześniej
stosowanymi w praktyce testami, takimi jak testy Kołmogorowa-Smirnowa,
Craméra von Misesa, [Liao i
Shimokawy](https://www.tandfonline.com/doi/abs/10.1080/00949659908811965)
oraz [Krita](http://www.numdam.org/item/JSFS_2014__155_3_135_0/).
Jednakże w mojej pracy dokonałam porównania tych testów na szerszym
zakresie alternatyw, które zostały podzielone według zachowania ich
funkcji hazardu.

Repozytorium zawiera nastepujące pliki:

- `newton_alg.R`: Implementuje algorytm Newtona-Raphsona do estymacji
  parametrów rozkładu Weibulla - $\lambda$ (parametr skali) i $\theta$
  (parametr kształtu).
- `test_form.R`: Generuje wartości wyżej wymienionych testów z
  pominięciem testu Krita.
- `Quantiles_generator_uncens.R`: Zawiera wszystkie niezbędne metody do
  generowania kwantyli dla pełnych prób.
- `powers_uncens.R`: Używa pętli `foreach` do generowania mocy testów na
  pełnej próbie, gdzie alternatywy były takie same jak w pracy
  wspomnianej powyżej, aby sprawdzić, czy uzyskane wyniki są podobne.
- `powers_uncens_HR.R` Używa pętli `foreach` do generowania mocy testów
  na pełnej próbie z alternatywami podzielonymi według zachowania ich
  funkcji hazardu.
- `powers_cens_HR.R`: Wykorzystuje metodę bootstrap do obliczania
  kwantyli dla prób cenzurowanych (z rozkładu jednostajnego) i pętle
  `foreach` do generowania mocy z alternatywami podzielonymi według
  zachowania ich funkcji hazardu.

Celem jest zbadanie efektywności nowych testów. W szczególności, pod
względem mocy, gdy ważnym czynnikiem jest typ funkcji hazardu.

Praca ta została napisana w języku polskim w 2024 r., kiedy byłam
studentką studiów pierwszego stopnia na Uniwersytecie Wrocławskim
(Instytut Matematyki i Informatyki).
