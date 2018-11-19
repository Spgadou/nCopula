#' Random number generator for Mother class objects
#'
#' @description Samples from a Mother class object.
#'
#' @param n the number of realisations.
#' @param structure an object of class Mother.
#'
#' @return A numeric matrix of sampled data from the structure
#'
#' @examples
#' ## Create the structure
#' structure <- GEO(0.1, 1, list(GAMMA(0.2, 2:3, NULL),
#'                         GEO(0.3, 4:5, NULL)))
#'
#' ## Sample from the structure
#' rCompCop(1000, structure)
#'
#' @importFrom stats rexp
#'
#' @author Simon-Pierre Gadoury
#' @export

rCompCop <- compiler::cmpfun(function(n, structure)
{
  e1 <- new.env()
  e1$res <- list()
  gen <- GeneticCodes(structure)
  e1$M0 <- structure@simul(n, as.numeric(structure@parameter))

  for (i in 1:length(gen))
  {
    if (length(gen[[i]]) == 2)
    {
      str2 <- Node(gen[[i]], structure)
      M.prec <- eval(parse(text = paste("e1$M", paste(gen[[i]][1], collapse = ""), sep = "")))

      R <- matrix(rexp(length(str2@arg) * n, 1), ncol = length(str2@arg), nrow = n)

      if (gen[[i]][length(gen[[i]])] == 0)
      {
        Theta <- M.prec

        ini <- stringr::str_replace_all(structure@Laplace, structure@Param, structure@parameter)
        ff <- function(z) eval(parse(text = ini))
        e1$res[[i]] <- ff(R / Theta)
      }
      else
      {
        Theta <- matrix(rep(vapply(1:length(M.prec), function(t) sum(str2@simul(M.prec[t], as.numeric(str2@parameter))), 0), length(str2@arg)),
                        ncol = length(str2@arg), nrow = n)

        ini <- stringr::str_replace_all(structure@PGF, structure@Param, structure@parameter)
        ini <- stringr::str_replace_all(ini, "z",
                                        stringr::str_replace_all(str2@Laplace, str2@Param, str2@parameter))
        ff <- function(z) eval(parse(text = ini))
        e1$res[[i]] <- ff(R / Theta)
      }
    }
    else if (length(gen[[i]]) > 2)
    {
      Lap <- stringr::str_replace_all(structure@PGF, structure@Param, structure@parameter)

      for (j in 2:(length(gen[[i]]) - 1))
      {
        str2 <- Node(gen[[i]][1:j], structure)

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

        variable0 <- paste("M", paste(gen[[i]][1:(j - 1)], collapse = ""), sep = "")
        variable1 <- paste("M", paste(gen[[i]][1:j], collapse = ""), sep = "")

        if (!exists(variable1, envir = e1))
        {
          eval(parse(text = paste("e1$", variable1, " <- vapply(1:length(", paste("e1$", variable0, sep = ""), "),
                                  function(t) sum(str2@simul(", paste("e1$", variable0, "[t],", sep = ""), "as.numeric(str2@parameter))), 0)", sep = "")))
        }
      }

      str2 <- Node(gen[[i]], structure)
      M.prec <- eval(parse(text = paste("e1$M", paste(gen[[i]][1:(length(gen[[i]]) - 1)], collapse = ""), sep = "")))

      if (gen[[i]][length(gen[[i]])] != 0)
      {
        ini <- stringr::str_replace_all(str2@Laplace, str2@Param, str2@parameter)
        Lap <- stringr::str_replace_all(Lap, "z", ini)

        Theta <- matrix(rep(vapply(1:length(M.prec), function(t) sum(str2@simul(M.prec[t], as.numeric(str2@parameter))), 0), length(str2@arg)),
                        ncol = length(str2@arg), nrow = n)
      }
      else
        Theta <- M.prec

      R <- matrix(rexp(length(str2@arg) * n, 1), ncol = length(str2@arg), nrow = n)

      ff <- function(z) eval(parse(text = Lap))
      e1$res[[i]] <- ff(R / Theta)
    }
  }
  do.call(cbind, e1$res)
})
