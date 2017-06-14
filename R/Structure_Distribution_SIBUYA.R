#' Construction of a SIBUYA Mother or Child Class Object
#'
#' @description Constructs either a Mother or Child class object for
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
#' @importFrom methods new
#' @examples
#' SIBUYA(0.5, NULL, list(GAMMA(1/30, c(5,6), NULL),
#'                        SIBUYA(0.1, NULL, list(GAMMA(1/30, c(1,2), NULL),
#'                                               GAMMA(1/30, c(3,4), NULL)))))
#' 
#' @author Simon-Pierre Gadoury
#' 
#' @export

SIBUYA <- compiler::cmpfun(function(par, unif, struc)
{
  if (length(unique(unif)) != length(unif))
    stop("The 'unif' argument must be composed of different values")

  if (is.null(struc))
  {
    t <- new("Sibuya_Child", parameter = par, arg = unif, dimension = length(unif), type = "Child", name = "Sibuya distribution", obj = "Sib")
  }

  else
  {
    if (class(struc) != "list")
      stop("The argument 'struc' must be a list")

    if (is.null(unif))
      t <- new("Sibuya_Mother", parameter = par, structure = struc, arg = 0, dimension = length(struc), type = "Mother", name = "Sibuya distribution", obj = "Sib")
    else
      t <- new("Sibuya_Mother", parameter = par, structure = struc, arg = unif, dimension = length(struc) + length(unif), type = "Mother", name = "Sibuya distribution", obj = "Sib")

  }

  t@Param <- "alpha"
  t@Laplace <- "1 - (1 - exp(-(z)))^(alpha)"
  t@LaplaceInv <- "-log(1 - (1 - (z))^(1/(alpha)))"
  t@PGF <- "1 - (1 - (z))^(alpha)"
  t@PGFInv <- "1 - (1 - (z))^(1/(alpha))"
  t@simul <- function(z, alpha) copula::rSibuya(z, alpha)
  t@theta <- vector("numeric")

  t
})
