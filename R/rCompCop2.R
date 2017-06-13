#' Density, Cdf, and Random Number Generator for Copulas Constructed Through Compounding
#'
#' @param n Number of realisations
#' @param str Object of class Mother
#'
#' @details rCompCop2 is more general (and easier to use) than rCompCop, but is slower.
#'
#' @author Simon-Pierre Gadoury
#'
#' @export

rCompCop2 <- compiler::cmpfun(function(n, str)
{
  e1 <- new.env()
  e1$res <- list()
  gen <- GeneticCodes(str)
  e1$M0 <- str@simul(n, str@parameter)

  for (i in 1:length(gen))
  {
    ## Ce bloc la est good je pense
    if (length(gen[[i]]) == 2)
    {
      str2 <- Node(gen[[i]], str)
      M.prec <- eval(parse(text = paste("e1$M", paste(gen[[i]][1], collapse = ""), sep = "")))
      Theta <- matrix(rep(vapply(1:length(M.prec), function(t) sum(str2@simul(M.prec[t], str2@parameter)), 0), length(str2@arg)),
                      ncol = length(str2@arg), nrow = n)
      R <- matrix(rexp(length(str2@arg) * n, 1), ncol = length(str2@arg), nrow = n)

      if (gen[[i]][length(gen[[i]])] == 0)
      {
        ini <- stringr::str_replace_all(str@Laplace, str@Param, str@parameter)
        ff <- function(z) eval(parse(text = ini))
        e1$res[[i]] <- ff(R / Theta)
      }
      else
      {
        ini <- stringr::str_replace_all(str@PGF, str@Param, str@parameter)
        ini <- stringr::str_replace_all(ini, "z",
                                        stringr::str_replace_all(str2@Laplace, str2@Param, str2@parameter))
        ff <- function(z) eval(parse(text = ini))
        e1$res[[i]] <- ff(R / Theta)
      }
    }
    else if (length(gen[[i]]) > 2)
    {
      ## Initialiser la Laplace du Theta final
      Lap <- stringr::str_replace_all(str@PGF, str@Param, str@parameter)

      for (j in 2:(length(gen[[i]]) - 1))
      {
        ## Sous structure
        str2 <- Node(gen[[i]][1:j], str)

        ## Conditions pour bâtir la Laplace
        if (gen[[i]][length(gen[[i]])] != 0)
        {
          ini <- stringr::str_replace_all(str2@PGF, str2@Param, str2@parameter)
          Lap <- stringr::str_replace_all(Lap, "z", ini)
        }
        else
        {
          if (j == length(gen[[i]]) - 1)
          {
            ini <- stringr::str_replace_all(str2@Laplace, str2@Param, str2@parameter)
            Lap <- stringr::str_replace_all(Lap, "z", ini)
          }
          else
          {
            ini <- stringr::str_replace_all(str2@PGF, str2@Param, str2@parameter)
            Lap <- stringr::str_replace_all(Lap, "z", ini)
          }
        }

        variable0 <- paste("M", paste(gen[[i]][1:(j - 1)], collapse = ""), sep = "") ## Le M au dessus du M
        variable1 <- paste("M", paste(gen[[i]][1:j], collapse = ""), sep = "") ## Le M

        if (!exists(variable1, envir = e1))
        {
          ## C'est important de bien le sentir
          eval(parse(text = paste("e1$", variable1, " <- vapply(1:length(", paste("e1$", variable0, sep = ""), "),
                                  function(t) sum(str2@simul(", paste("e1$", variable0, "[t],", sep = ""), "str2@parameter)), 0)", sep = "")))
        }
      }

      str2 <- Node(gen[[i]], str)
      M.prec <- eval(parse(text = paste("e1$M", paste(gen[[i]][1:(length(gen[[i]]) - 1)], collapse = ""), sep = "")))

      if (gen[[i]][length(gen[[i]])] != 0)
      {
        ini <- stringr::str_replace_all(str2@Laplace, str2@Param, str2@parameter)
        Lap <- stringr::str_replace_all(Lap, "z", ini)
      }

      Theta <- matrix(rep(vapply(1:length(M.prec), function(t) sum(str2@simul(M.prec[t], str2@parameter)), 0), length(str2@arg)),
                      ncol = length(str2@arg), nrow = n)
      R <- matrix(rexp(length(str2@arg) * n, 1), ncol = length(str2@arg), nrow = n)

      ff <- function(z) eval(parse(text = Lap))
      e1$res[[i]] <- ff(R / Theta)
    }
  }
  do.call(cbind, e1$res)
})
