#' Density, Cdf, and Random Number Generator for Copulas
#'
#' @param n Number of realisations
#' @param copula An Archimedean copula class object
#' @return Sampled data.
#'
#' @author Simon-Pierre Gadoury
#'
#' @export

rCop <- compiler::cmpfun(function(n, copula)
{
  sim <- copula@theta
  param <- copula@parameter
  dim <- copula@dimension

  param2 <- copula@depend
  param2 <- stringr::str_replace_all(param2, "alpha", param)
  param2 <- eval(parse(text = param2))

  phi <- stringr::str_replace_all(copula@phi, "alpha", param)
  phi <- parse(text = stringr::str_replace_all(phi, "z", "res"))

  res <- t(-log(runif(dim*n)) / matrix(sim(n, param2), ncol = n, nrow = dim, byrow = TRUE))

  eval(phi)
})
