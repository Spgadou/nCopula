#' Cdf Function for Archimedean Copulas (archm) Class Objects
#'
#' @description The cdf function for Archimedean copulas (arcm) class objetcs.
#' @param copula An Archimedean copula (archm) class object
#' @param vector If false, returns a function with (u_1, u_2, ..., u_dim) as arguments, else,
#' just (u)
#' @return A function ...
#'
#' @author Simon-Pierre Gadoury
#'
#' @examples
#' cop <- Clayton(5, 2)
#' pCop(cop, vector = TRUE)(c(0.1, 1))
#' pCop(cop, vector = FALSE)(0.1, 0.1)
#'
#' @seealso \link{rCop}, \link{Clayton}, \link{AMH}, \link{Gumbel}, \link{Frank}
#'
#' @export

pCop <- compiler::cmpfun(function(copula, vector = TRUE, express = FALSE)
{
  phi <- copula@phi
  dim <- copula@dimension
  phi.inv <- copula@phi.inv

  if (vector)
    uu <- paste("u[", 1:dim, "]", sep = "")
  else
    uu <- paste("u", 1:dim, sep = "")

  res <- numeric(dim)
  for (i in 1:dim)
    res[i] <- stringr::str_replace_all(phi.inv, "z", uu[i])
  res <- paste("(", res, ")", collapse = " + ")
  cop <- stringr::str_replace_all(phi, "z", res)

  t1 <- "function(z)"

  if (vector)
    t3 <- paste(c("u"), collapse = ", ")
  else
  {
    tt <- paste(uu, collapse = ", ")
    t2 <- paste(c(tt), collapse = ", ")
    t3 <- t2
  }

  expr2 <- "eval(parse(text = cop))"

  input <- stringr::str_replace_all(t1, "z", t3)
  input2 <- paste(c(input, expr2), collapse = " ")

  alpha <- copula@parameter

  res2 <- parse(text = input2)

  if (express == FALSE)
    return(eval(res2))
  else
    return(cop)
})
