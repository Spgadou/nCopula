#' Construction of a GAMMA Child Class Object
#'
#' @description The function GAMMA constructs a gamma Child class object for
#' a given parameter and arguments.
#'
#' @param par parameter of the distribution.
#' @param unif uniform structure, a numeric vector of grouped
#' numbers, i.e. c(1,2,3) is translated as being c(u1, u2, u3).
#' @param structure nesting structure of the form
#'
#' X(par1, c(i,...), list(Y(par2, c(j,...), NULL),
#'                        Z(par3, c(k,...), NULL))),
#'
#' where X, Y, and Z are compatible functions (see 'details').
#' It is to note that if structure is NULL, the function will automatically
#' be of class Child. For continuous distributions (i.e. GAMMA), structure is
#' always NULL.
#'
#' @author Simon-Pierre Gadoury
#'
#' @family mother or child class objects.
#'
#' @importFrom methods new
#'
#' @examples
#' GEO(0.5, NULL, list(GAMMA(1/30, c(5,6), NULL),
#'                     GEO(0.1, NULL, list(GAMMA(1/30, c(1,2), NULL),
#'                                         GAMMA(1/30, c(3,4), NULL)))))
#' @export

GAMMA <- compiler::cmpfun(function(par, unif, structure = NULL)
{
     if (length(unique(unif)) != length(unif))
          stop("The 'unif' argument must be composed of different values")

     if (par < 0)
          stop("Wrong 'param' input")

     if (!is.null(structure))
          stop("Argument 'structure' must be NULL for a 'Child' class")

     t <- new("Gamma_Child", parameter = as.character(par), arg = unif, type = "Child", dimension = length(unif), name = "Gamma distribution", obj = "Gamma")

     t@Param <- "alpha"
     t@Laplace <- "(1 / (1 + (z)))^(alpha)"
     t@LaplaceInv <- "((z)^(-1/(alpha)) - 1)"
     t@simul <- function(z, alpha) rgamma(z, alpha, 1)
     t@theta <- vector("numeric")
     t@Der <- function(tt, k, type)
     {
          if (type == "Laplace")
          {
            if (k == 0)
              stringr::str_replace_all(t@Laplace, "z", tt)
            else
            {
             expr1 <- paste("(", 0:(k - 1), " + alpha)", collapse = " * ", sep = "")
             ini <- paste("(-1)^(k) * ", expr1, " * (1 + (z))^(-alpha - (k))", sep = "")
             ini <- stringr::str_replace_all(ini, "z", tt)
             stringr::str_replace_all(ini, "k", as.character(k))
            }
          }
          else if (type == "LaplaceInv")
          {
               stringr::str_replace_all("-1/(alpha) * (z)^(-1/(alpha)-1)", "z", tt)
          }
     }
     t@FUN <- function(type)
     {
          if (type == "Laplace")
               function(tt, alpha) (1 + (tt))^(-alpha)
          else if (type == "LaplaceInv")
               function(tt, alpha) (tt)^(-1/alpha) - 1
          else if (type == "Laplace.Der")
          {
               function(tt, alpha, k, expon = 1)
               {
                    if (expon == 0)
                         0
                    else
                    {
                         (-1)^k * prod((0:(k-1) + (alpha * expon))) * (1 + (tt))^(-(alpha * expon) - k)
                    }
               }
          }
          else if (type == "LaplaceInv.Der")
          {
               function(tt, alpha)
                    -(1/alpha) * (tt)^(-1/alpha - 1)
          }
     }

     t
})
