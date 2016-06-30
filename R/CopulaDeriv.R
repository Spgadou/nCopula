#' Copula derivative
#'
#' Take the derivative of any copula.
#' @param nest The copula function
#' @param parameters Number of copulas (if nested)
#' @return The symbolical derivative
#' @export

CopulaDeriv <- function(nest, parameters = 1)
{
  ui <- head(names(formals(nest)), -parameters)
  der <- list()
  der[[1]] <- Deriv(nest, ui[1], cache.exp = FALSE)
  for (i in 2:length(ui))
    der[[i]] <- Deriv(der[[i - 1]], ui[i], cache.exp = FALSE)
  res <- Simplify(der[[length(ui)]])
  res
}
