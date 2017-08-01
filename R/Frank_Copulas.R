#' Construction of an Archimedean Copula Class Object
#'
#' Constructs Frank Archimedean copula object with a given parameter and dimension.
#'
#' @description Constructs a Frank Archimedean copula object with
#' a given parameter and dimension.
#' @param dim Dimension of the copula (>= 2), which is, by default, 2.
#' @param param Parameter of the copula.
#' @param density Compute the expression of the density of the copulas.
#'
#' @author Simon-Pierre Gadoury
#'
#' @return An archm S4 class object.
#'
#' @importFrom copula rlog
#' @export

Frank <- compiler::cmpfun(function(param, dim = 2L)
{
  if (param < 0)
    stop("Wrong 'param' input")

  verif <- eval(parse(text = paste0(dim, "L")))

  if (!is.integer(verif) || dim <= 1)
    stop("The dimension must be an integer greater than or equal to 2")

  phi <- "log(1 - (alpha) * exp(-(z))) / log(1 - (alpha))"
  phi.inv <- "-log((1 - (1 - (alpha))^(z)) / (alpha))"
  dep.param <- "alpha"
  rBiv <- function(n, alpha, u) -log(1 - (runif(n) * (exp(-alpha) - 1)) / (runif(n) * (exp(-alpha * u) - 1) - exp(-alpha * u))) / alpha
  th <- function(z, alpha) copula::rlog(z, alpha)

  new("frank",
      phi = phi,
      phi.inv = phi.inv,
      theta = th,
      rBiv = rBiv,
      depend = dep.param,
      dimension = dim,
      parameter = param,
      name = "Frank copula")
  if (density)
  {
       tt <- GAMMA(1/10, 1:dim, NULL)
       
       uu <- paste("u", 1:dim, sep = "")
       expr1 <- numeric(dim)
       for (i in 1:dim)
            expr1[i] <- stringr::str_replace_all(tt@Der("z", 1, "LaplaceInv"), "z", uu[i])
       expr1 <- paste("(", expr1, ")", sep = "", collapse = " * ")
       
       nu <- numeric(dim)
       for(i in 1:dim)
            nu[i] <- stringr::str_replace_all(tt@LaplaceInv, "z", uu[i])
       nu <- paste("(", nu, ")", sep = "", collapse = " + ")
       
       expr2 <- stringr::str_replace_all(tt@Der("z", dim, "Laplace"), "z", nu)
       densit <- paste("(", expr1, ") * (", expr2, ")", sep = "")
       densit <- stringr::str_replace_all(densit, "alpha", "(1/alpha)")
       
       new("frank",
           phi = phi,
           phi.inv = phi.inv,
           rBiv = rBiv,
           theta = th,
           depend = dep.param,
           dimension = dim,
           parameter = param,
           dens = densit,
           name = "Frank copula")
  }
  else
  {
       new("frank",
           phi = phi,
           phi.inv = phi.inv,
           rBiv = rBiv,
           theta = th,
           depend = dep.param,
           dimension = dim,
           parameter = param,
           name = "Frank copula")
  }
})
