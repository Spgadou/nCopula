#' Construction of a Mother or Child Class Object
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
#' @author Simon-Pierre Gadoury
#'
#' @importFrom methods new
#' @examples
#' GEO(0.5, NULL, list(GAMMA(1/30, c(5,6), NULL),
#'                     GEO(0.1, NULL, list(GAMMA(1/30, c(1,2), NULL),
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

#' Construction of a Mother or Child Class Object
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
        ini <- stringr::str_replace_all("factorial(k) / (uu)^(k - 1) / gamma * ((z) / (uu))^2 * ((z)/((uu) * gamma) - 1)^(k - 1)", "z",
                                        t@PGF)
        ini <- stringr::str_replace_all(ini, "k", k)
        ini <- stringr::str_replace_all(ini, "uu", tt)
        stringr::str_replace_all(ini, "z", tt)
      }
      else if (type == "PGFInv")
      {
        ini <- stringr::str_replace_all("(1 - gamma) * (z)^2", "z",
                                        t@PGFInv)
        stringr::str_replace_all(ini, "z", tt)
      }
    }
  }
  t@simul <- function(z, gamma) rgeom(z, gamma) + 1
  t@theta <- vector("numeric")
  t@cop <- function(gamma, dim) AMH(gamma, dim)
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

  t
})

#' Construction of a Child Class Object
#'
#' @description Constructs a Child class object for
#' a given parameter and arguments
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
#' @importFrom methods new
#' @examples
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

#' Construction of a Mother or Child Class Object
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
#' @author Simon-Pierre Gadoury
#'
#' @importFrom methods new
#' @examples
#' GEO(0.5, NULL, list(GAMMA(1/30, c(5,6), NULL),
#'                     GEO(0.1, NULL, list(GAMMA(1/30, c(1,2), NULL),
#'                                         GAMMA(1/30, c(3,4), NULL)))))
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
