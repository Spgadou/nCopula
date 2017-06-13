#' Construction of the copula with a known structure
#'
#' @param str Object of class Mother
#' @param lvl Used internally
#' @param j used internally
#'
#' @details This is a test.
#'
#' @author Simon-Pierre Gadoury
#'
#' @export

rCompCop2 <- function(str, lvl = 0, j = 1)
{
  if (str@type == "Mother")
  {
    if (lvl == 0)
    {
      gen <- GeneticCodes(str)
      e1 <- new.env()
      e1$C <- stringr::str_replace_all(str@PGF, str@Param, str@parameter)
      e1$M0 <- numeric(str@dimension)
    }
    else
    {
      ini <- stringr::str_replace_all(str@PGF, str@Param, str@parameter)
      eval(parse(text = paste("e1$M", lvl - 1, "[j] <- ini", sep = "")))
      eval(parse(text = paste("e1$M", lvl, " <- numeric(length(str@dimension))", sep = "")))
    }

    for (i in 1:str@dimension)
    {
      FUN(str@structure[[i]], lvl + 1, i)
    }
  }
  else
  {
    ini <- stringr::str_replace_all(str@Laplace, str@Param, str@parameter)
    eval(parse(text = paste("e1$M", lvl - 1, "[j] <- ini", sep = "")))
  }
}





