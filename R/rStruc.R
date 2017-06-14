#' Sampling for compound rvs
#'
#' @description The function rStruc() generate n samples from the compound rvs.
#' @param n Number of realisations
#' @param str Object of class Mother (the structure)
#'
#' @example rStruc(10000, GEO(0.5, NULL, list(GAMMA(1/30, c(5,6), NULL),
#'                                      GEO(0.1, NULL, list(GAMMA(1/30, c(1,2), NULL),
#'                                      GAMMA(1/30, c(3,4), NULL))))))
#'
#' @return Return...
#' @author Simon-Pierre Gadoury
#'
#' @export

rStruc <- compiler::cmpfun(function(n, str)
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

      if (gen[[i]][length(gen[[i]])] == 0)
      {
        Theta <- M.prec
        e1$res[[i]] <- Theta
      }
      else
      {
        Theta <- matrix(rep(vapply(1:length(M.prec), function(t) sum(str2@simul(M.prec[t], str2@parameter)), 0), length(str2@arg)),
                        ncol = length(str2@arg), nrow = n)
        e1$res[[i]] <- Theta
      }
    }
    else if (length(gen[[i]]) > 2)
    {
      for (j in 2:(length(gen[[i]]) - 1))
      {
        ## Sous structure
        str2 <- Node(gen[[i]][1:j], str)

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
        Theta <- matrix(rep(vapply(1:length(M.prec), function(t) sum(str2@simul(M.prec[t], str2@parameter)), 0), length(str2@arg)),
                        ncol = length(str2@arg), nrow = n)
      }
      else
        Theta <- M.prec

      e1$res[[i]] <- Theta
      }
  }
  do.call(cbind, e1$res)
})
