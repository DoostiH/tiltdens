# tiltdens

<!-- badges: start -->
[![R-CMD-check](https://github.com/DoostiH/tiltdens/actions/workflows/R-CMD-check.yaml/badge.svg)](https://github.com/DoostiH/tiltdens/actions/workflows/R-CMD-check.yaml)
[![License: MIT](https://img.shields.io/badge/license-MIT-blue.svg)](LICENSE.md)
<!-- badges: end -->

Tilted and data-sharpened nonparametric density estimation in R.

A conventional kernel density estimator with a non-negative kernel cannot
converge faster than $O_p(n^{-4/5})$, however smooth the underlying density is.
Higher-order kernels beat that rate, but the price is an estimate that goes
negative and oscillates in the tails — it is no longer a density, and it
suggests structure that is not there.

This package implements a way out. Perturb a conventional kernel estimator so
that it sits as close as possible to a higher-order one, either by re-weighting
the observations (*tilting*) or by moving them (*data sharpening*). Because the
kernel is a proper density and the weights are a probability vector, the result
is always non-negative and always integrates to one, while inheriting the fast
convergence rate of the estimator it was fitted against.

The methods are from:

- Doosti, H. and Hall, P. (2016). Making a non-parametric density estimator more
  attractive, and more accurate, by data perturbation.
  *Journal of the Royal Statistical Society B* **78**, 445–462.
  [doi:10.1111/rssb.12112](https://doi.org/10.1111/rssb.12112)
- Doosti, H., Hall, P. and Mateu, J. (2018). Nonparametric tilted density
  function estimation: a cross-validation criterion.
  *Journal of Statistical Planning and Inference* **197**, 51–68.
  [doi:10.1016/j.jspi.2017.12.003](https://doi.org/10.1016/j.jspi.2017.12.003)

## Installation

```r
# install.packages("remotes")
remotes::install_github("DoostiH/tiltdens", build_vignettes = TRUE)
```

The only dependency beyond base R is **quadprog**.

## Usage

```r
library(tiltdens)

set.seed(2016)
x <- c(rnorm(50, -1.5), rnorm(50, 1.5))

fit <- tilt_density_cv(x)     # 2018 method: fast, no comparator needed
fit
plot(fit)
```

Fitted objects inherit from `"density"`, so `plot()`, `lines()` and `points()`
work exactly as they do for `stats::density()`.

The point of the exercise, in two numbers:

```r
min(sinc_density(x)$y)   # the infinite-order comparator dips below zero
#> [1] -0.0162

min(fit$y)               # a tilted estimate cannot
#> [1] 5.1e-07
```

## Which function do I want?

| | Function | Use it when |
|---|---|---|
| **Start here** | `tilt_density_cv()` | You want a good estimate quickly. Chooses bandwidth and weights together by cross-validation; no comparator needed. |
| | `tilt_density()` | You want the estimate anchored to a specific infinite-order estimator, or to compare the sinc and trapezoidal comparators. |
| | `sharpen_density()` | You want the observations moved rather than re-weighted, and the shifts are of interest. Slowest of the three. |
| | `sinc_density()`, `trapezoid_density()` | You want the infinite-order comparators on their own. |

Bandwidths come from `bw_comparator_cv()` (the default), `bw_flattop()` and
`bw_nrd_robust()`.

## Kernels

The papers use the standard normal throughout, noting that it "simplified
numerical work". The theory in Section 3.1 of the 2016 paper is actually stated
for kernels satisfying its condition (3.1), which holds when `K` is a *k*-fold
convolution of a Laplace density, so that family is available too.

```r
tilt_kernels()
#> [1] "gaussian"     "laplace"      "laplace2"     "laplace3"     "laplace4"
#> [6] "epanechnikov" "biweight"     "triweight"    "triangular"

tilt_density(x, m = 3, kernel = "laplace2")   # the paper's k = 2 example
```

Custom kernels can be supplied as a list of `dens`, `ft` and `support`.

Weight selection stays a convex quadratic program whatever the kernel: the Gram
matrix is `A[i,j] = (K * K)(x_i - x_j)`, whose Fourier transform is `phi^2 >= 0`,
so `A` is positive semidefinite by construction. Equivalently, `t(p) %*% A %*% p`
is `integral of fhat^2`, which cannot be negative.

Bandwidths mean different things for different kernels — `bw = 0.5` smooths far
less with the Epanechnikov, supported on `[-1,1]`, than with the Gaussian. The
package handles this with the canonical factor of Marron and Nolan (1988),
`(R(K)/mu2(K)^2)^(1/5)`, so the default bandwidth is rescaled to whichever
kernel you pick:

```r
bw_canonical_factor("epanechnikov")          # 1.7188
bw_convert(0.5, from = "gaussian", to = "epanechnikov")
```

Changing the kernel therefore changes the *shape* of the fit rather than how
much it is smoothed. On a bimodal sample, the spread in integrated squared
error across four kernels falls from 2.3x without the rescaling to 1.2x with it.
For the Gaussian the factor is one, so the papers' setting is unchanged.

## Multivariate data

Pass a matrix. Section 3.3 of the 2016 paper replaces the kernel with a product
kernel and keeps a scalar bandwidth, which makes every quantity the method needs
factorise across dimensions — so the l-variate matrices are elementwise products
of the univariate ones, and the optimisation problem is unchanged.

```r
X <- cbind(c(rnorm(60, -1.5), rnorm(60, 1.5)),
           c(rnorm(60, -1),   rnorm(60, 1)))
fit <- tilt_density(X)
plot(fit)                       # contour, with the observations overlaid
predict(fit, cbind(0, 0))
```

The same thing happens in two dimensions as in one: the sinc comparator goes
negative over parts of the plane, and the tilted estimate cannot.

Weights are per-observation in the multivariate case (`m` does not apply — block
structure rests on the order statistics, which have no natural counterpart above
one dimension). Data sharpening stays univariate: the search would range over
`n * l` variables, which is impractical rather than merely slow.

## Shape constraints

The estimator is *linear* in the weights, so a shape requirement evaluated on a
grid is a set of linear inequalities, and imposing it leaves the problem a
convex quadratic program.

```r
x <- c(rnorm(40), rnorm(6, 2.6, 0.12))    # a spurious secondary bump
tilt_density(x, constraint = "unimodal")  # also "increasing", "decreasing"
```

Shape-constrained tilting is due to Hall and Huang (2002), which both papers
cite. Applying the same device to the 2018 cross-validation criterion follows
the extension suggested in that paper's discussion.

Constraints need enough degrees of freedom to work with: with only two or three
distinct weights there may be no feasible solution, and the error message says
so. Use `m = Inf`.

`vignette("tiltdens")` walks through the whole thing.

## How many distinct weights?

The `m` argument controls how many distinct values the weights may take.
`m = Inf` gives every observation its own weight; `m = 3` allows three, over a
central block and two tails.

More freedom is not automatically better. In the papers' simulations `m = 3` was
often *more* accurate, because fewer free parameters means less overfitting to
the comparator's own noise.

When `m` is finite, *where* the blocks begin and end matters too. The `breaks`
argument offers three strategies:

```r
c(equal   = tilt_density(x, m = 3, breaks = "equal")$distance2,
  modal   = tilt_density(x, m = 3, breaks = "modal")$distance2,
  optimal = tilt_density(x, m = 3, breaks = "optimal")$distance2)
#>    equal    modal  optimal
#> 0.001975 0.002135 0.001629
```

`"optimal"` is the default for `tilt_density()` and treats the boundaries as
part of the optimisation, as Section 4.1 of the 2016 paper specifies. On bimodal
data it finds the trough between the modes on its own. `"equal"` is Algorithm A
of the 2018 paper and is much cheaper. `"modal"` places boundaries at the
troughs of a pilot estimate. You can also pass the boundaries yourself.

## When these methods help

They earn their keep on complex densities: sharp peaks, well-separated modes,
heavy tails. On a simple smooth density a conventional kernel estimator is hard
to beat, and these methods will roughly match it rather than improve on it. The
papers' simulations show exactly that pattern, and it is worth knowing before
you reach for them.

## Notes on the implementation

Both papers' criteria reduce to the **same convex quadratic program** over the
probability simplex, differing only in one linear term: the 2016 method compares
against a comparator estimator, the 2018 method against the leave-one-out fit.
One solver therefore serves both, and the solution is unique rather than
something a search has to hunt for.

A few things differ from the original MATLAB code written for the papers:

- Weight selection is solved exactly with `quadprog` (with a projected-gradient
  fallback), rather than by a derivative-free pattern search over a
  100-dimensional simplex.
- The quadrature for the oscillatory cross terms adapts its node count to the
  frequency. A fixed rule loses all accuracy once the pairwise data distance
  exceeds a few tens of bandwidths.
- Block boundaries are computed in one place from the ranks of the sample, so
  the grouping and the weight assignment cannot drift apart.
- Comparator bandwidths are computed from the data rather than loaded from
  precomputed files.

`NEWS.md` has the full list.

The test suite checks mathematical identities against independent computations —
that the Gram matrix really equals the squared L2 norm of the estimator, that
the cross terms match numerical integration for both comparator kernels, that
the cross-validation criterion matches its own definition evaluated by brute
force — rather than checking that output matches stored numbers.

## MATLAB

The original MATLAB code from the papers, cleaned up and documented, lives at
[DoostiH/tilted-density-estimation](https://github.com/DoostiH/tilted-density-estimation).
The two implementations agree to 12 significant figures on every shared
quantity.

## Citation

```r
citation("tiltdens")
```

Please cite both papers as well as the software.

## Contributing

See [CONTRIBUTING.md](CONTRIBUTING.md). Two extensions would be particularly
welcome: the multivariate case of Section 3.3 of the 2016 paper, and shape
constraints imposed through a generalised cross-validation criterion, as
suggested in the discussion of the 2018 paper.

## License

MIT © Hassan Doosti
