#' Construction of an Archimedean Copula Class Object.
#'
#' Constructs an AMH Archimedean copula object with a given parameter and dimension.
#'
#' @description Constructs an AMH Archimedean copula object with
#' a given parameter and dimension.
#' @param dim dimension of the copula (>= 2), which is, by default, 2.
#' @param param parameter of the copula.
#' @param density compute the expression of the density of the copulas.
#'
#' @author Simon-Pierre Gadoury
#'
#' @importFrom stats rgeom
#'
#' @return An archm S4 class object.
#'
#' @export

AMH <- compiler::cmpfun(function(param, dim = 2L, density = FALSE)
{
     if (param < 0 || param >= 1)
          stop("Wrong 'param' input")

     verif <- eval(parse(text = paste0(dim, "L")))

     if (!is.integer(verif) || dim <= 1)
          stop("The dimension must be an integer greater than or equal to 2")

     phi <- "(alpha) / (exp(z) - (1 - alpha))"
     phi.inv <- "log((alpha + (z) * (1 - alpha)) / (z))"
     dep.param <- "alpha"
     param.th <- "alpha"
     th <- function(z, alpha) rgeom(z, alpha) + 1

     param <- as.character(param)

     if (density)
     {
          tt <- GEO(1/10, 1:dim, NULL)

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
          densit <- stringr::str_replace_all(densit, "gamma", "alpha")

          new("amh",
              phi = phi,
              phi.inv = phi.inv,
              theta = th,
              depend = dep.param,
              dimension = dim,
              parameter = param,
              dens = densit,
              par.th = param.th,
              name = "AMH copula")
     }
     else
     {
          new("amh",
              phi = phi,
              phi.inv = phi.inv,
              theta = th,
              depend = dep.param,
              dimension = dim,
              parameter = param,
              par.th = param.th,
              name = "AMH copula")
     }
})
