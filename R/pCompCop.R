#' Construction of the copula with a known structure
#'
#' @param str Object of class Mother
#'
#' @return An expression in terms of u
#'
#' @author Simon-Pierre Gadoury
#'
#' @export

pCompCop <- function(str)
{
  e1 <- new.env(hash = TRUE, parent = parent.frame(), size = 10L)
  e1$gen <- GeneticCodes(str)
  str_ini <- str

  FUN <- function(str, lvl = 0, j = 1, v = 0)
  {
    if (str@type == "Mother")
    {
      if (lvl == 0)
      {
        e1$C <- stringr::str_replace_all(str@PGF, str@Param, str@parameter)
        e1$M0 <- numeric(str@dimension - length(str@arg) + 1 * (sum(str@arg) != 0))
      }
      else
      {
        ini <- stringr::str_replace_all(str@PGF, str@Param, str@parameter)
        eval(parse(text = paste("e1$M", lvl - 1, "[j] <- ini", sep = "")))
        eval(parse(text = paste("e1$M", lvl, " <- numeric(str@dimension - length(str@arg) + length(str@arg) * (sum(str@arg) != 0))", sep = "")))
      }

      for (i in 1:(str@dimension - length(str@arg)))
      {
        FUN(str@structure[[i]], lvl + 1, i, v = c(v, i))
      }

      if (sum(str@arg) != 0)
      {
        charr <- InvLap_Child(c(v, 0), str_ini)
        uu <- paste("u", str@arg, sep = "")
        res <- numeric(length(uu))
        for (i in 1:length(uu))
          res[i] <- stringr::str_replace_all(charr, "z", uu[i])
        res <- paste("(", res, ")", collapse = " * ")

        eval(parse(text = paste("e1$M", lvl, "[str@dimension - length(str@arg) + 1] <- res", sep = "")))
      }

      char1 <- paste("(", eval(parse(text = paste("e1$M", lvl, sep = ""))), ")", collapse = " * ")

      if (lvl > 0)
      {
        eval(parse(text = paste("e1$M", lvl - 1, "[j] <- stringr::str_replace_all(", paste("e1$M", lvl - 1, "[j]", sep = ""), ", 'z', '", char1, "')", sep = "")))
      }
      else
        e1$C <- stringr::str_replace_all(e1$C, "z", char1)
    }
    else
    {
      argum <- str@arg
      uu <- paste("u", argum, sep = "")
      nu <- InvLap_Child(v, str_ini)
      res <- numeric(length(argum))
      for (y in 1:length(argum))
        res[y] <- stringr::str_replace_all(nu, "z", uu[y])
      res <- paste("(", res, ")", collapse = " + ")

      ini <- stringr::str_replace_all(str@Laplace, str@Param, str@parameter)
      ini <- stringr::str_replace_all(ini, "z", res)
      eval(parse(text = paste("e1$M", lvl - 1, "[j] <- ini", sep = "")))
    }
  }

  FUN(str)
  e1$C
}








