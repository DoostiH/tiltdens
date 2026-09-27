## Test environments
* local Windows, R 4.3.3
* win-builder, R-devel
* GitHub Actions: Windows, macOS, Ubuntu (release, devel, oldrel-1)

## R CMD check results
0 errors | 0 warnings | 0 notes

## This release
Version 0.2.0 responds to an editorial review of the accompanying Journal of
Statistical Software submission. It fixes a bug in bw_flattop() for data with
a large scale (affecting multivariate input such as datasets::faithful), adds a
bounded and pluggable optimizer to sharpen_density(), adds make_kernel() and
check_kernel(), adds conventional_density(), and corrects documentation. The
sharpen_density() example now runs in about two seconds. See NEWS.md.
