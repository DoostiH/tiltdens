# JSS manuscript

`tiltdens-jss.Rmd` is a skeleton for a Journal of Statistical Software
submission. It has the section structure JSS expects and, in HTML comments, what
belongs in each section and why. Replace the comments with prose.

## Building it

```r
install.packages("rticles")
rmarkdown::render("tiltdens-jss.Rmd")
```

`rticles::jss_article` supplies `jss.cls` and the required formatting. JSS also
accepts `.Rnw`; if you prefer Sweave, the same structure carries over.

## House style

JSS is strict about markup, and reviewers do send papers back over it:

- Software in `\pkg{}`, functions in `\code{}`, languages in `\proglang{}`.
- Code chunks use the `R> ` prompt, which the setup chunk already configures.
- Everything must be reproducible from the sources you submit, including the
  simulation. Set seeds.
- The paper and the package are reviewed together, so keep them in step: if a
  function is renamed, the manuscript changes too.

## The argument to lead with

The two papers are already published, so a restatement will not carry a software
paper. The contribution to foreground is that both criteria reduce to the *same*
convex quadratic program over the probability simplex, differing in one linear
term. That single observation is what makes the package possible: it gives a
unique solution instead of a derivative-free search over a 100-dimensional
simplex, and it is what lets block boundaries, shape constraints and the
multivariate case be added without changing the solver.

Section 3 is where that lives. Give it room.

## One thing to be ready for

A reviewer may try to reproduce Table 1 of either paper and find the numbers
differ. Get ahead of it. The reasons are known and documented in the MATLAB
repository: the original bandwidths were selected from among several local
minima of a multimodal cross-validation criterion in a way that is not
recoverable, and a defect in the original sampler affected two of the eight test
densities. A short paragraph saying so, in Section 5 or the discussion, is far
better than being asked.
