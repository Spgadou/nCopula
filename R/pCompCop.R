#' Distribution function of Mother class objects
#'
#' @description Distribution function of a Mother class object.
#'
#' @param structure object of class Mother.
#' @param vector logical. If false, returns a function or a character string with (u_1, u_2, ...) as arguments, else,
#' just (u).
#' @param express logical. If false, returns a function, else, a character string.
#'
#' @return The distribution function in the form of either a function or a character string.
#'
#' @examples
#' ## Create the structure
#' structure <- LOG(0.5, NULL, list(GAMMA(1/30, c(5,6), NULL),
#'                               LOG(0.1, NULL, list(GAMMA(1/30, c(1,2), NULL),
#'                               GAMMA(1/30, c(3,4), NULL)))))
#'
#' ## Character string
#' pCompCop(structure, vector = TRUE, express = TRUE)
#' pCompCop(structure, vector = FALSE, express = TRUE)
#'
#' ## Function
#' pCompCop(structure, vector = TRUE, express = FALSE)
#' pCompCop(structure, vector = FALSE, express = FALSE)
#'
#' @export

pCompCop <- function(structure, vector = FALSE, express = TRUE)
{
  e1 <- new.env(hash = TRUE, parent = parent.frame(), size = 10L)
  e1$gen <- GeneticCodes(structure)
  e1$argmax <- 0
  str_ini <- structure

  FUN <- function(structure, lvl = 0, j = 1, v = 0)
  {
    if (structure@type == "Mother")
    {
      argum <- structure@arg
      for (kk in 1:length(argum))
      {
        if (argum[kk] > e1$argmax)
          e1$argmax <- argum[kk]
      }

      if (lvl == 0)
      {
        e1$C <- stringr::str_replace_all(structure@PGF, structure@Param, structure@parameter)
        e1$M0 <- numeric(structure@dimension - length(structure@arg) + 1 * (sum(structure@arg) != 0))
      }
      else
      {
        ini <- stringr::str_replace_all(structure@PGF, structure@Param, structure@parameter)
        eval(parse(text = paste("e1$M", lvl - 1, "[j] <- ini", sep = "")))
        eval(parse(text = paste("e1$M", lvl, " <- numeric(structure@dimension - length(structure@arg) + length(structure@arg) * (sum(structure@arg) != 0))", sep = "")))
      }

      for (i in 1:(structure@dimension - structure@arg))
      {
        FUN(structure@structure[[i]], lvl + 1, i, v = c(v, i))
      }

      if (sum(structure@arg) != 0)
      {
        charr <- InvLap(c(v, 0), str_ini)
        uu <- paste("u", structure@arg, sep = "")
        res <- numeric(length(uu))
        for (i in 1:length(uu))
          res[i] <- stringr::str_replace_all(charr, "z", uu[i])
        res <- paste("(", res, ")", collapse = " * ")

        eval(parse(text = paste("e1$M", lvl, "[structure@dimension - length(structure@arg) + 1] <- res", sep = "")))
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
      argum <- structure@arg
      for (kk in 1:length(argum))
      {
        if (argum[kk] > e1$argmax)
          e1$argmax <- argum[kk]
      }

      uu <- paste("u", argum, sep = "")
      nu <- InvLap(v, str_ini)
      res <- numeric(length(argum))
      for (y in 1:length(argum))
        res[y] <- stringr::str_replace_all(nu, "z", uu[y])
      res <- paste("(", res, ")", collapse = " + ")

      ini <- stringr::str_replace_all(structure@Laplace, structure@Param, structure@parameter)
      ini <- stringr::str_replace_all(ini, "z", res)
      eval(parse(text = paste("e1$M", lvl - 1, "[j] <- ini", sep = "")))
    }
  }

  FUN(structure)

  cop <- e1$C
  dim <- e1$argmax

  if (express)
  {
    if (vector)
    {
      for (i in dim:1)
        cop <- stringr::str_replace_all(cop, paste("u", i, sep = ""), paste("u[", i, "]", sep = ""))
      cop
    }
    else
      cop
  }
  else
  {
    if (vector)
    {
      for (i in dim:1)
        cop <- stringr::str_replace_all(cop, paste("u", i, sep = ""), paste("u[", i, "]", sep = ""))
      eval(parse(text = "function(u) eval(parse(text = cop))"))
    }
    else
    {
      ff <- "function(z) eval(parse(text = cop))"
      uu <- paste("u", 1:dim, sep = "", collapse = ", ")
      eval(parse(text = stringr::str_replace_all(ff, "z", uu)))
    }
  }
}




