# JSS manuscript

| File | What it is |
|---|---|
| `tiltdens-jss.pdf` | The manuscript, 12 pages, compiled |
| `tiltdens-jss.Rnw` | Its source (Sweave + `jss.cls`) |
| `tiltdens.bib` | References |
| `code.R` | Standalone replication script — reproduces every result, table and figure |
| `sim/` | Cached simulation results, read by the `.Rnw` when building the PDF |
| `jss.cls`, `jss.bst`, `jsslogo.jpg` | JSS style files, so the PDF builds anywhere |

## What JSS asks for at submission

From their author information page:

1. **PDF manuscript in JSS style** — `tiltdens-jss.pdf`. Done.
2. **Source code for the software** — the package. On CRAN once accepted; the
   tarball otherwise.
3. **Replication materials** — `code.R`, plus a file `code.html` produced by
   running `knitr::spin("code.R")`, which must end with `sessionInfo()`. The
   script already ends that way. **Generate `code.html` on your own machine**,
   because JSS wants the session information from a real platform, and because
   it takes about fifteen minutes:

   ```r
   install.packages("knitr")   # if needed
   setwd("path/to/paper")
   knitr::spin("code.R")
   ```

   This runs every line and writes `code.html`. Compare its numbers against the
   PDF; they should match, with at most last-digit differences in the
   optimised-boundary fits.

Their other stated requirements, and where we stand:

- *Discuss advantages and disadvantages against existing implementations, with
  empirical illustrations* — Sections 5 and 6, and Table 2.
- *Extensive simulation studies are discouraged* — Section 5.1 is four densities
  and forty replicates, and says why.
- *Replication within one hour on a regular PC* — about fifteen minutes.
- *S3 classes with print, plot and summary methods* — all three exist, plus
  `predict`. The `summary` method was added for this submission.
- *GPL-compatible licence* — MIT is GPL-compatible.
- *Software on CRAN, not just GitHub* — submitted 9 September 2026, pending
  manual inspection.
- *They encourage including the JSS article as a package vignette* — optional;
  worth doing after acceptance, when the text is final.

## Rebuilding the PDF

```r
Sweave("tiltdens-jss.Rnw")
```

then `pdflatex`, `bibtex`, `pdflatex`, `pdflatex` on `tiltdens-jss.tex`. The
simulation table is read from `sim/`, so this takes under a minute. To
regenerate `sim/` from scratch, run `code.R`; its Section 5.1 block produces the
same numbers.

## Two things to check before submitting

**The existing-software paragraph in Section 1.2.** It claims no package in any
language implements tilted density estimation. That was true when written;
check it still is.

**Table 2.** The capability table for `ks`, `KernSmooth`, `kdensity`,
`logcondens` and `sharpData` was compiled from their documentation. Worth a
quick re-read of each package's current help before submitting.
