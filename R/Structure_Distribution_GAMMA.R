#' Construction of a GAMMA Child Class Object
#'
#' @description The function GAMMA constructs a gamma Child class object for
#' a given parameter and arguments.
#'
#' @param par Dimension of the distribution
#' @param unif Uniform structure, a numeric vector of grouped
#' numbers.
#'
#' i.e. c(1,2,3) is translated as being c(u1, u2, u3).
#' @param struc Nesting structure of the form
#'
#' X(par1, c(i,...), list(Y(par2, c(j,...), NULL),
#'                        Z(par3, c(k,...), NULL))),
#'
#' where X, Y, and Z are compatible functions (see 'details').
#' It is to note that if struc is NULL, the function will automatically
#' be of class Child. For continuous distributions (i.e. GAMMA), struc is
#' always NULL.
#'
#' @slot Param The name of the parameter used
#' @slot parameter The value of the parameter
#' @slot dimension The dimension
#' @slot type The type of function (either child or mother)
#' @slot arguments The corresponding arguments (ex.: arguments 1 and 2 imply 'u1' and 'u2')
#' @slot structure The structure below the node of type 'Mother'
#' @slot Laplace Expression of the LST
#' @slot LaplaceInv Expression of the inverse LST
#' @slot PGF Expression of the pgf
#' @slot PGFInv Expression of the inverse pgf
#' @slot simul Fonction to sample from the distribution
#' @slot theta I don't know honestly
#' @slot cop Construct an Archimedean copula with this distribution
#' @slot Der Fonction to compute the expression of the 'k'th derivative of either the 'PGF', 'PGFInv', 'Laplace' or 'LaplaceInv'
#' @slot FUN Fonction to compute the function of the 'k'th derivative of either the 'PGF', 'PGFInv', 'Laplace' or 'LaplaceInv'
#'
#' @author Simon-Pierre Gadoury
#'
#' @family mother or child class objects
#'
#' @importFrom methods new
#'
#' @examples
#' We recall that since GAMMA is continuous, the structure is exclusively a child.
#' GEO(0.5, NULL, list(GAMMA(1/30, c(5,6), NULL),
#'                     GEO(0.1, NULL, list(GAMMA(1/30, c(1,2), NULL),
#'                                         GAMMA(1/30, c(3,4), NULL)))))
#' @export

GAMMA <- compiler::cmpfun(function(par, unif, struc = NULL)
{
     if (length(unique(unif)) != length(unif))
          stop("The 'unif' argument must be composed of different values")

     if (par < 0)
          stop("Wrong 'param' input")

     if (!is.null(struc))
          stop("Argument 'struc' must be NULL for a 'Child' class")

     t <- new("Gamma_Child", parameter = par, arg = unif, type = "Child", dimension = length(unif), name = "Gamma distribution", obj = "Gamma")

     t@Param <- "alpha"
     t@Laplace <- "(1 / (1 + (z)))^(alpha)"
     t@LaplaceInv <- "((z)^(-1/(alpha)) - 1)"
     t@simul <- function(z, alpha) rgamma(z, alpha, 1)
     t@theta <- vector("numeric")
     t@Der <- function(tt, k, type)
     {
          if (type == "Laplace")
          {
               expr1 <- paste("(", 0:(k - 1), " + alpha)", collapse = " * ", sep = "")
               ini <- paste("(-1)^(k) * ", expr1, " * (1 + (z))^(-alpha - (k))", sep = "")
               ini <- stringr::str_replace_all(ini, "z", tt)
               stringr::str_replace_all(ini, "k", k)
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
