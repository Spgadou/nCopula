#' Construction of an Archimedean Copula Class Object
#'
#' Constructs an AMH Archimedean copula object with a given parameter and dimension
#'
#' @param dim Dimension of the copula (>= 2), which is, by default, 2
#' @param param Parameter of the copula
#' @param density Should the expression of the density be computed ?
#'
#' @author Simon-Pierre Gadoury
#'
#' @return An archm S4 class object
#'
#' @export

AMH <- compiler::cmpfun(function(param, dim = 2L)
{
     if (param < 0 || param >= 1)
          stop("Wrong 'param' input")
     
     verif <- eval(parse(text = paste0(dim, "L")))
     
     if (!is.integer(verif) || dim <= 1)
          stop("The dimension must be an integer greater than or equal to 2")
     
     phi <- "(alpha) / (exp(z) - (1 - alpha))"
     phi.inv <- "log((alpha + (z) * (1 - alpha)) / (z))"
     dep.param <- "alpha"
     th <- function(z, alpha) rgeom(z, alpha) + 1
     
     new("amh",
         phi = phi,
         phi.inv = phi.inv,
         theta = th,
         depend = dep.param,
         dimension = dim,
         parameter = param,
         name = "AMH copula")
})