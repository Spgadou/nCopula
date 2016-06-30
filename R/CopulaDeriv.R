#' Copula derivative
#'
#' Take the derivative of any copula.
#' @param cop The copula function
#' @param parameters Number of copulas (if nested)
#' @return The symbolical derivative
#' @export

CopulaDeriv <- function(cop, parameters = 1)
{
  ui <- head(names(formals(nest)), -parameters)
  der <- list()
  der[[1]] <- Deriv(nest, ui[1], cache.exp = FALSE)
  for (i in 2:length(ui))
    der[[i]] <- Deriv(der[[i - 1]], ui[i], cache.exp = FALSE)
  res <- Simplify(der[[length(ui)]])
  res
}

#' Bivariate Clayton Copula
#'
#' Bivariate clayton copula function
#' @param alpha Dependence parameter
#' @return The clayton copula
#' @export

BiClayton <- function(u1, u2, alpha) (u1^(-alpha) + u2^(-alpha) - 1)^(-1/alpha)

#' Bivariate Copula Derivative
#'
#' Take the derivative of bivariate copulas.
#' @param cop The copula function
#' @return The symbolical derivative
#' @details Clayton, Frank.
#' @export

BivDeriv <- function(cop)
{
  Deriv(Deriv(cop, "u1", cache.exp = FALSE), "u2", cache.exp = FALSE)
}
