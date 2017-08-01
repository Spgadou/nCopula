#' Construction of an Archimedean Copula Class Object
#'
#' Constructs Frank Archimedean copula object with a given parameter and dimension
#'
#' @param dim Dimension of the copula (>= 2), which is, by default, 2
#' @param param Parameter of the copula
#' @param density Should the expression of the density be computed ?
#'
#' @author Simon-Pierre Gadoury
#'
#' @return An archm S4 class object
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
})
