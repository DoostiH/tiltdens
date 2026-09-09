## Test environments
* local Windows, R 4.3.3
* win-builder, R-devel (2026-09-08 r90509 ucrt)
* GitHub Actions: Windows, macOS, Ubuntu (release, devel, oldrel-1)

## R CMD check results
0 errors | 0 warnings | 1 note

* This is a new release.
* The words flagged as possibly misspelled in DESCRIPTION are the surnames
  Doosti and Mateu, and the technical terms "comparator" and "sinc" (the sinc
  kernel, whose Fourier transform is flat on [-1, 1]). All are correct as
  written.

## Notes
The package implements two published methods (Doosti and Hall 2016, JRSS-B;
Doosti, Hall and Mateu 2018, JSPI) that had no implementation available in any
language. Both DOIs are cited in the DESCRIPTION.