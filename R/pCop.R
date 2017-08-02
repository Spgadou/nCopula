#' Distribution function of archm class objects
#'
#' @description Distribution function of an Archimedean copula (archm) class object.
#'
#' @param copula an Archimedean copula (archm) class object.
#' @param vector logical. If false, returns a function or a character string with (u_1, u_2, ..., u_dim) as arguments, else,
#' just (u).
#' @param express logical. If false, returns a function, else, a character string.
#'
#' @return The distribution function in the form of either a function or a character string.
#'
#' @author Simon-Pierre Gadoury
#'
#' @examples
#' cop <- Clayton(5, 2)
#' pCop(cop, vector = TRUE, express = TRUE)
#' pCop(cop, vector = FALSE, express = TRUE)
#'
#' @seealso \link{rCop}, \link{Clayton}, \link{AMH}, \link{Gumbel}, \link{Frank}
#'
#' @export

pCop <- compiler::cmpfun(function(copula, vector = FALSE, express = TRUE)
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
    eval(res2)
  else
    cop
})
