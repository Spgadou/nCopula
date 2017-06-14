#' Cdf, and Random Number Generator for Copulas
#'
#' @param copula An Archimedean copula class object
#' @param vector If false, returns a function with (u_1, u_2, ..., u_dim) as arguments, else, just (u)
#' @return Either an expression, function, code, or sampled data.
#'
#' @author Simon-Pierre Gadoury
#'
#' @export

pCop <- compiler::cmpfun(function(copula, vector = TRUE)
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
  }

  expr2 <- "eval(parse(text = cop))"

  input <- stringr::str_replace_all(t1, "z", t3)
  input2 <- paste(c(input, expr2), collapse = " ")

  alpha <- copula@parameter

  res2 <- parse(text = input2)
  return(eval(res2))
})
