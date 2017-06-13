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

pCompCop2 <- function(str)
{
  e1 <- new.env(hash = TRUE, parent = parent.frame(), size = 10L)
  e1$gen <- GeneticCodes(str)

  FUN <- function(str, lvl = 0, j = 1)
  {
    if (str@type == "Mother")
    {
      if (lvl == 0)
      {
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

      char1 <- paste("(", eval(parse(text = paste("e1$M", lvl, sep = ""))), ")", collapse = " * ")

      if (lvl > 0)
        eval(parse(text = paste("e1$M", lvl - 1, "[j] <- stringr::str_replace_all(", paste("e1$M", lvl - 1, "[j]", sep = ""), ", 'z', ", char1, ")", sep = "")))
      else
        e1$C <- stringr::str_replace_all(e1$C, "z", char1)
    }
    else
    {
      ini <- stringr::str_replace_all(str@Laplace, str@Param, str@parameter)
      eval(parse(text = paste("e1$M", lvl - 1, "[j] <- ini", sep = "")))
    }
  }

  FUN(str)
  e1$C
}

pCompCop2(str2)





