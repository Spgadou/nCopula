#' Construction of a GEO Mother or Child Class Object
#'
#' @description Constructs either a GEO Mother or Child class object for
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
#' @family mother or child class objects
#'
#' @author Simon-Pierre Gadoury
#'
#' @importFrom methods new
#' @examples
#' GEO(0.5, NULL, list(GAMMA(1/30, c(5,6), NULL),
#'                     GEO(0.1, NULL, list(GAMMA(1/30, c(1,2), NULL),
#'                                         GAMMA(1/30, c(3,4), NULL)))))
#' @export

GEO <- compiler::cmpfun(function(par, unif, struc)
{
     if (length(unique(unif)) != length(unif))
          stop("The 'unif' argument must be composed of different values")

     if (par > 1 || par < 0)
          stop("Wrong 'param' input")

     if (is.null(struc))
     {
          t <- new("Geo_Child", parameter = par, arg = unif, dimension = length(unif), type = "Child", name = "Shifted geometric distribution", obj = "Geo")
     }

     else
     {
          if (class(struc) != "list")
               stop("The argument 'struc' must be a list")

          if (is.null(unif))
               t <- new("Geo_Mother", parameter = par, structure = struc, arg = 0, dimension = length(struc), type = "Mother", name = "Shifted geometric distribution", obj = "Geo")
          else
               t <- new("Geo_Mother", parameter = par, structure = struc, arg = unif, dimension = length(struc) + length(unif), type = "Mother", name = "Shifted geometric distribution", obj = "Geo")

     }

     if (t@type == "Mother")
     {
          t@Param <- "gamma"
          t@Laplace <- "(gamma)*exp(-(z)) / (1 - (1 - (gamma)) * exp(-(z)))"
          t@LaplaceInv <- "-log(1 / (((gamma)/(z)) + (1 - (gamma))))"
          t@PGF <- "(gamma)*(z) / (1 - (1-(gamma))*(z))"
          t@PGFInv <- "1 / (((gamma)/(z)) + (1 - (gamma)))"
          t@Der <- function(tt, k, type)
          {
               if (type == "PGF")
               {
                    if (k >= 1)
                    {
                         ini <- stringr::str_replace_all("factorial(k) / (uu)^(k - 1) / gamma * ((z) / (uu))^2 * ((z)/((uu) * gamma) - 1)^(k - 1)", "z",
                                                         t@PGF)
                         ini <- stringr::str_replace_all(ini, "k", k)
                         ini <- stringr::str_replace_all(ini, "uu", tt)
                         stringr::str_replace_all(ini, "z", tt)
                    }
                    else
                         t@PGF
               }
               else if (type == "PGFInv")
               {
                    ini <- "(gamma) / (gamma + (z) * (1 - gamma))^2"
                    stringr::str_replace_all(ini, "z", tt)
               }
          }
          t@FUN <- function(type)
          {
               if (type == "PGF")
                    function(tt, gamma) gamma * (tt) / (1 - (1 - gamma) * (tt))
               else if (type == "PGFInv")
                    function(tt, gamma) (tt) / (gamma + (tt) * (1 - gamma))
               else if (type == "PGF.Der")
               {
                    function(tt, gamma, k)
                         factorial(k) / (tt)^(k - 1) / gamma * (t@FUN("PGF")(tt, gamma) / (tt))^2 * (t@FUN("PGF")(tt, gamma) / ((tt) * gamma) - 1)^(k - 1)
               }
               else if (type == "PGFInv.Der")
               {
                    function(tt, gamma)
                         gamma / (gamma + (tt) * (1 - gamma))^2
               }
          }
     }
     else
     {
          t@Param <- "alpha"
          t@Laplace <- "(alpha)*exp(-(z)) / (1 - (1 - (alpha)) * exp(-(z)))"
          t@LaplaceInv <- "-log(1 / (((alpha)/(z)) + (1 - (alpha))))"
          t@PGF <- "(alpha)*(z) / (1 - (1-(alpha))*(z))"
          t@PGFInv <- "1 / (((alpha)/(z)) + (1 - (alpha)))"
          t@Der <- function(tt, k, type)
          {
               if (type == "PGF")
               {
                    ini <- stringr::str_replace_all("factorial(k) / (uu)^(k - 1) / alpha * ((z) / (uu))^2 * ((z)/((uu) * alpha) - 1)^(k - 1)", "z",
                                                    t@PGF)
                    ini <- stringr::str_replace_all(ini, "k", k)
                    ini <- stringr::str_replace_all(ini, "uu", tt)
                    stringr::str_replace_all(ini, "z", tt)
               }
               else if (type == "PGFInv")
               {
                    ini <- stringr::str_replace_all("(1 - alpha) * (z)^2", "z",
                                                    t@PGFInv)
                    stringr::str_replace_all(ini, "z", tt)
               }
          }
     }
     t@simul <- function(z, gamma) rgeom(z, gamma) + 1
     t@theta <- vector("numeric")
     t@cop <- function(gamma, dim) AMH(gamma, dim)

     t
})
