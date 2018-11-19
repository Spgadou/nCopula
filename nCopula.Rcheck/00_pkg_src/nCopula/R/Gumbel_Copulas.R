#' Construction of an Archimedean Copula Class Object
#'
#' @description Constructs a Gumbel Archimedean copula object with
#' a given parameter and dimension.
#'
#' @param dim dimension of the copula (>= 2), which is, by default, 2
#' @param param parameter of the copula
#'
#' @author Simon-Pierre Gadoury
#'
#' @return An archm S4 class object.
#'
#' @importFrom copula rstable1
#'
#' @export

Gumbel <- compiler::cmpfun(function(param, dim = 2L)
{
     param.th <- NULL

     if (param < 1)
          stop("Wrong 'param' input")

     verif <- eval(parse(text = paste0(dim, "L")))

     if (!is.integer(verif) || dim <= 1)
          stop("The dimension must be an integer greater than or equal to 2")

     phi <- "exp(-(z)^(1/alpha))"
     phi.inv <- "(-log(z))^alpha"
     dep.param <- "alpha"
     par.th <- "alpha"
     th <- function(z, alpha) copula::rstable1(z, 1/alpha, 1, cos(pi/(2*alpha))^alpha, 0, 1)

     param <- as.character(param)

     new("gumbel",
         phi = phi,
         phi.inv = phi.inv,
         theta = th,
         depend = dep.param,
         dimension = dim,
         parameter = param,
         par.th = param.th,
         name = "Gumbel copula")
})
