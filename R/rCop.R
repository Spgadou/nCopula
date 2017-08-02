#' Random number generator for Archimedean copula class objects
#'
#' @description Random number generator for archm class objects.
#'
#' @param n number of realisations.
#' @param copula an Archimedean copula (archm) class object.
#'
#' @details For bivariate archm copula objects, the function uses the conditional approach.
#' As for dimensions higher than 2, the Marshall-Olkin (1988) approach is chosen instead.
#'
#' @return A numeric matrix containing the samples.
#'
#' @examples
#' ## Create the trivariate archm copula object
#' cop <- Clayton(5, 3)
#'
#' ## Generate the samples
#' res <- rCop(10000, cop)
#'
#' ## Plot the values
#' pairs(res, pch = 16, cex = 0.7)
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

  param2 <- copula@par.th
  param2 <- stringr::str_replace_all(param2, "alpha", param)
  param2 <- eval(parse(text = param2))

  phi <- stringr::str_replace_all(copula@phi, "alpha", as.character(param))
  phi <- parse(text = stringr::str_replace_all(phi, "z", "res"))

  res <- t(-log(runif(dim*n)) / matrix(sim(n, param2), ncol = n, nrow = dim, byrow = TRUE))

  eval(phi)
})
