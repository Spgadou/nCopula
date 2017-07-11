#' Random number generator for Archimedean copula class objects
#'
#' @description The function rCop is a random number generator for archm class object.
#'
#' @param n Number of realisations
#' @param copula An Archimedean copula (archm) class object
#'
#' @return A matrix containing the samples
#'
#' @examples
#' cop <- Clayton(5, 3)
#'
#' rCop(1000, cop)
#'
#' @seealso \link{pCop}, \link{Clayton}, \link{AMH}, \link{Frank}, \link{Gumbel}
#'
#' @author Simon-Pierre Gadoury
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
