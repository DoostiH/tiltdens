# Contributing

Bug reports, questions and pull requests are all welcome.

## Reporting a problem

Open an issue with a small reproducible example — ideally a few lines that
someone else can paste into a fresh R session. Include the output of
`sessionInfo()`.

If a density estimate looks wrong, it helps a great deal to say what you
expected and why. Some surprising output is correct: the comparator estimators
are *supposed* to go negative, and the cross-validation criterion is *supposed*
to be multimodal.

## Pull requests

* Add a test for anything you change. The suite lives in `tests/testthat/` and
  runs with `devtools::test()`.
* Tests here check mathematical identities against independent computations
  rather than against stored output. If you add a new formula, please check it
  against numerical integration or a closed form, as the existing tests do.
* Run `devtools::check()` before opening the request.
* Document exported functions with roxygen2 and regenerate with
  `devtools::document()`.

## Scope

The package covers univariate tilting and data sharpening. Two extensions would
be welcome and are not yet implemented:

* The multivariate case described in Section 3.3 of Doosti and Hall (2016).
* Shape constraints such as unimodality or bimodality, imposed through a
  generalised cross-validation criterion, as suggested in the discussion of
  Doosti, Hall and Mateu (2018).
