#' Construction of a LOG Mother or Child Class Object
#'
#' @description Constructs either a LOG Mother or Child class object for
#' a given parameter, arguments, and nesting structure.
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
#' @family mother or child class objects
#'
#' @author Simon-Pierre Gadoury
#'
#' @importFrom methods new
#' @examples
#'LOG(0.5, NULL, list(GAMMA(1/30, c(5,6), NULL), 
#'                     LOG(0.1, NULL, list(GAMMA(1/30, c(1.2), NULL), 
#'                                         GAMMA(1/30, c(3,4), NULL)))))
#' @export

LOG <- compiler::cmpfun(function(par, unif, struc)
{
     if (length(unique(unif)) != length(unif))
          stop("The 'unif' argument must be composed of different values")
     
     if (par > 1 || par < 0)
          stop("Wrong 'param' input")
     
     if (is.null(struc))
     {
          t <- new("Log_Child", parameter = par, arg = unif, dimension = length(unif), name = "Logarithmic distribution", type = "Child", obj = "Log")
     }
     
     else
     {
          if (class(struc) != "list")
               stop("The argument 'struc' must be a list")
          
          if (is.null(unif))
               t <- new("Log_Mother", parameter = par, dimension = length(struc), structure = struc, arg = 0, name = "Logarithmic distribution", type = "Mother", obj = "Log")
          else
               t <- new("Log_Mother", parameter = par, dimension = length(struc) + length(unif), structure = struc, arg = unif, name = "Logarithmic distribution", type = "Mother", obj = "Log")
     }
     
     if (t@type == "Mother")
     {
          t@Param <- "gamma"
          t@Laplace <- "log(1 - (gamma) * exp(-(z))) / log(1 - (gamma))"
          t@LaplaceInv <- "-log((1 - (1 - (gamma))^(z)) / (gamma))"
          t@PGF <- "log(1 - (gamma)*(z)) / log(1 - (gamma))"
          t@PGFInv <- "((1 - (1 - (gamma))^(z))/(gamma))"
     }
     else
     {
          t@Param <- "alpha"
          t@Laplace <- "log(1 - (alpha) * exp(-(z))) / log(1 - (alpha))"
          t@LaplaceInv <- "-log((1 - (1 - (alpha))^(z)) / (alpha))"
          t@PGF <- "log(1 - (alpha)*(z)) / log(1 - (alpha))"
          t@PGFInv <- "((1 - (1 - (alpha))^(z))/(alpha))"
     }
     t@simul <- function(z, gamma) copula::rlog(z, gamma)
     t@theta <- vector("numeric")
     t@cop <- function(gamma, dim) Frank(gamma, dim)
     
     t
})