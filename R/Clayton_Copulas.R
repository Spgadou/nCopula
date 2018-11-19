#' Construction of an Archimedean Copula Class Object
#'
#' @description Constructs a Clayton Archimedean copula object with
#' a given parameter and dimension.
#'
#' @param dim the dimension of the copula (>= 2), which is, by default, 2.
#' @param param the parameter of the copula.
#' @param density logical. Should the expression of the density of the copula be computed?
#'
#' @return An archm S4 class object.
#'
#' @importFrom stats rgamma runif
#'
#' @author Simon-Pierre Gadoury
#'
#' @export

Clayton <- compiler::cmpfun(function(param, dim = 2L, density = FALSE)
{
     if (param < 0)
          stop("Wrong 'param' input")

     verif <- eval(parse(text = paste0(dim, "L")))

     if (!is.integer(verif) || dim <= 1)
          stop("The dimension must be an integer greater than or equal to 2")

     phi <- "exp(log((z) + 1)*(-1/alpha))"
     phi.inv <- "((z)^(-alpha) - 1)"
     dep.param <- "alpha"
     param.th <- "(1/alpha)"
     rBiv <- function(n, alpha, u) (u^(-alpha) * (runif(n)^(-alpha / (alpha + 1)) - 1) + 1)^(-1/alpha)
     th <- function(z, alpha) rgamma(z, alpha, 1)

     param <- as.character(param)

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
          densit <- stringr::str_replace_all(densit, "alpha", dep.param)

          new("clayton",
              phi = phi,
              phi.inv = phi.inv,
              rBiv = rBiv,
              theta = th,
              depend = dep.param,
              dimension = dim,
              parameter = param,
              dens = densit,
              par.th = param.th,
              name = "Clayton copula")
     }
     else
     {
          new("clayton",
              phi = phi,
              phi.inv = phi.inv,
              rBiv = rBiv,
              theta = th,
              depend = dep.param,
              dimension = dim,
              parameter = param,
              par.th = param.th,
              name = "Clayton copula")
     }
})
