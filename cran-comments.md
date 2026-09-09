## Test environments

* local: Ubuntu 24.04, R 4.3.3
* (add win-builder and R-hub results before submitting)

## R CMD check results

0 errors | 0 warnings | 1 note

* This is a new release.

## Notes for the reviewer

The package implements two published estimators (Doosti and Hall 2016,
JRSS-B; Doosti, Hall and Mateu 2018, JSPI) that had no implementation in any
language available to users. Both DOIs are cited in the DESCRIPTION.

Examples and tests are kept small so they run quickly. `sharpen_density()` is
the only slow function; its example uses a deliberately small optimisation
budget, which is documented in the help page.
