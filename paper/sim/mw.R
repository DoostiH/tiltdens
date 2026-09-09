## Marron & Wand (1992) normal-mixture densities 1-8
mw <- list(
  Gau = list(w = 1,               mu = 0,                  sd = 1),
  SkU = list(w = c(1,1,3)/5,      mu = c(0, 1/2, 13/12),   sd = c(1, 2/3, 5/9)),
  StS = list(w = rep(1/8, 8),     mu = 3*((2/3)^(0:7) - 1), sd = (2/3)^(0:7)),
  KtU = list(w = c(2,1)/3,        mu = c(0, 0),            sd = c(1, 1/10)),
  Out = list(w = c(1,9)/10,       mu = c(0, 0),            sd = c(1, 1/10)),
  Bim = list(w = c(1,1)/2,        mu = c(-1, 1),           sd = c(2/3, 2/3)),
  SeB = list(w = c(1,1)/2,        mu = c(-3/2, 3/2),       sd = c(1/2, 1/2)),
  SkB = list(w = c(3,1)/4,        mu = c(0, 3/2),          sd = c(1, 1/3))
)
mw_pdf <- function(d) function(t) Reduce(`+`, Map(function(w, m, s) w * dnorm(t, m, s), d$w, d$mu, d$sd))
mw_rnd <- function(d, n) { k <- sample(length(d$w), n, TRUE, d$w); rnorm(n, d$mu[k], d$sd[k]) }
