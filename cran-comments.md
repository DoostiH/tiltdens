## Test environments
* local Windows, R 4.6.1
* win-builder, R-devel
* GitHub Actions: Windows, macOS, Ubuntu (release, devel, oldrel-1)

## R CMD check results
0 errors | 0 warnings | 1 note

* Days since last update: this is a small update requested in the editorial
  review of the accompanying Journal of Statistical Software submission. It
  adds validation of user-supplied kernels, so that an unsuitable kernel can
  no longer be used for estimation by mistake, and summary() and plot()
  methods for kernel objects. See NEWS.md.
