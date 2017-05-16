#' Bivariate Dependent Sum Approximation (CDF)
#'
#' Bivariate lower and upper bound approximation
#' @param s CDF evaluation
#' @param m Precision
#' @param FUN The copula
#' @param dep Dependence parameter
#' @return The lower and upper bound approximation of the CDF at s
#' @export

repBiv <- function(s, m, FUN, dep)
{
  options(error = function() stop("Problem between the chair and the screen"))

  rect.lower <- 2^m - 1
  rect.upper <- 2^m

  lower <- numeric(rect.lower)
  for (i in 1:rect.lower)
    lower[i] <- FUN(c(i/(2^m) * s, (2^m - i)/(2^m) * s), dep) - FUN(c((i - 1)/(2^m) * s, (2^m - i)/(2^m) * s), dep)

  upper <- numeric(rect.upper)
  for (i in 1:rect.upper)
    upper[i] <- FUN(c(i/(2^m) * s, (2^m - i + 1) /(2^m) * s), dep) - FUN(c((i - 1)/(2^m) * s, (2^m - i + 1) /(2^m) * s), dep)

  l <- sum(lower)
  u <- sum(upper)

  if (l > 1)
    l <- 1
  if (u > 1)
    u <- 1

  list("Lower bound" = l, "Upper bound" = u)
}
