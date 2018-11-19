#' Construction of a LOG Mother or Child Class Object
#'
#' @description Constructs either a LOG Mother or Child class object for
#' a given parameter, arguments, and nesting structure.
#'
#' @param par parameter of the distribution.
#' @param unif uniform structure, a numeric vector of grouped
#' numbers, i.e. c(1,2,3) is translated as being c(u1, u2, u3).
#' @param structure nesting structure of the form
#'
#' X(par1, c(i,...), list(Y(par2, c(j,...), NULL),
#'                        Z(par3, c(k,...), NULL))),
#'
#' where X, Y, and Z are compatible functions (see 'details').
#' It is good to note that if structure is NULL, the function will automatically
#' be of class Child. For continuous distributions (i.e. GAMMA), structure is
#' always NULL.
#'
#' @family mother or child class objects.
#'
#' @author Simon-Pierre Gadoury
#'
#' @importFrom methods new
#'
#' @examples
#' LOG(0.5, NULL, list(GAMMA(1/30, c(5,6), NULL),
#'                     LOG(0.1, NULL, list(GAMMA(1/30, c(1.2), NULL),
#'                                         GAMMA(1/30, c(3,4), NULL)))))
#' @export

LOG <- compiler::cmpfun(function(par, unif, structure)
{
     tt <- NULL

     if (length(unique(unif)) != length(unif))
          stop("The 'unif' argument must be composed of different values")

     if (par > 1 || par < 0)
          stop("Wrong 'param' input")

     if (is.null(structure))
     {
          t <- new("Log_Child", parameter = as.character(par), arg = unif, dimension = length(unif), name = "Logarithmic distribution", type = "Child", obj = "Log")
     }

     else
     {
          if (class(structure) != "list")
               stop("The argument 'structure' must be a list")

          if (is.null(unif))
               t <- new("Log_Mother", parameter = as.character(par), dimension = length(structure), structure = structure, arg = 0, name = "Logarithmic distribution", type = "Mother", obj = "Log")
          else
               t <- new("Log_Mother", parameter = as.character(par), dimension = length(structure) + length(unif), structure = structure, arg = unif, name = "Logarithmic distribution", type = "Mother", obj = "Log")
     }

     if (t@type == "Mother")
     {
          t@Param <- "gamma"
          t@Laplace <- "log(1 - (gamma) * exp(-(z))) / log(1 - (gamma))"
          t@LaplaceInv <- "-log((1 - (1 - (gamma))^(z)) / (gamma))"
          t@PGF <- "log(1 - (gamma)*(z)) / log(1 - (gamma))"
          t@PGFInv <- "((1 - (1 - (gamma))^(z))/(gamma))"
          t@Der <- function(tt, k, type)
          {
            if (type == "PGF")
            {
              if (k >= 1)
              {
                ini <- paste("(", -factorial(k - 1), ") / log(1 - gamma) * (gamma / (1 - gamma * (z)))^(", k, ")", sep = "")
                stringr::str_replace_all(ini, "z", tt)
              }
              else
                t@PGF
            }
            else if (type == "PGFInv")
            {
              if (k == 1)
              {
                ini <- "-log(1 - gamma) / gamma * (1 - gamma)^(z)"
                stringr::str_replace_all(ini, "z", tt)
              }
              else if (k == 0)
                t@PGFInv
            }
          }
          t@FUN <- function(type)
          {
            if (type == "PGF")
              function(tt, gamma) log(1 - (gamma)*(tt)) / log(1 - (gamma))
            else if (type == "PGFInv")
              function(tt, gamma) ((1 - (1 - (gamma))^(tt))/(gamma))
            else if (type == "PGF.Der")
            {
              function(tt, gamma, k)
                -factorial(k - 1) / log(1 - gamma) * (gamma / (1 - gamma * (tt)))^(k)
            }
            else if (type == "PGFInv.Der")
            {
              function(tt, gamma)
                -log(1 - gamma) / gamma * (1 - gamma)^(tt)
            }
          }
     }
     else
     {
          t@Param <- "alpha"
          t@Laplace <- "log(1 - (alpha) * exp(-(z))) / log(1 - (alpha))"
          t@LaplaceInv <- "-log((1 - (1 - (alpha))^(z)) / (alpha))"
          t@PGF <- "log(1 - (alpha)*(z)) / log(1 - (alpha))"
          t@PGFInv <- "((1 - (1 - (alpha))^(z))/(alpha))"
          t@Der <- function(tt, k, type)
          {
            if (type == "PGF")
            {
              if (k >= 1)
              {
                ini <- paste("(", -factorial(k - 1), ") / log(1 - alpha) * (alpha / (1 - alpha * (z)))^(", k, ")", sep = "")
                stringr::str_replace_all(ini, "z", tt)
              }
              else
                stringr::str_replace_all(t@PGF, "z", tt)
            }
            else if (type == "PGFInv")
            {
              if (k == 1)
              {
                ini <- "-log(1 - alpha) / alpha * (1 - alpha)^(z)"
                stringr::str_replace_all(ini, "z", tt)
              }
              else if (k == 0)
                stringr::str_replace_all(t@PGFInv, "z", tt)
            }
            else if (type == "Laplace")
            {
              if (k >= 1)
              {
                res <- numeric(k)
                for (r in 1:k)
                {
                  ini <- t@Der("exp(-(z))", r, "PGF")
                  ini <- paste("(", ini, ") * (exp(-", r, " * (z)) * (-1)^(", k, "))", sep = "")

                  input <- (-1)^(r - 1) / factorial(1) / factorial(r - 1) * 1^k
                  if (r > 1)
                  {
                    for (s in 2:r)
                    {
                      input <- input + ((-1)^(r - s) / factorial(s) / factorial(r - s) * s^k)
                    }
                  }

                  res[r] <- paste("(", ini, ") * (", input, ")", sep = "")
                }
                res <- paste("(", res, ")", sep = "", collapse = " + ")
                stringr::str_replace_all(res, "z", tt)
              }
              else
              {
                stringr::str_replace_all(t@Laplace, "z", tt)
              }
            }
            else if (type == "LaplaceInv")
            {
              ini <- paste("-(", t@Der(tt, 1, "PGFInv"), ") / (", stringr::str_replace_all(t@PGFInv,
                                                                                           "z",
                                                                                           tt), ")", sep = "")
              ini
            }
          }
          t@FUN <- function(type)
          {
            if (type == "PGF")
              function(tt, gamma) log(1 - (gamma)*(tt)) / log(1 - (gamma))
            else if (type == "PGFInv")
              function(tt, gamma) ((1 - (1 - (gamma))^(tt))/(gamma))
            else if (type == "Laplace")
              log(1 - (gamma) * exp(-(tt))) / log(1 - (gamma))
            else if (type == "LaplaceInv")
              -log((1 - (1 - (gamma))^(tt)) / (gamma))
            else if (type == "PGF.Der")
            {
              function(tt, gamma, k)
                -factorial(k - 1) / log(1 - gamma) * (gamma / (1 - gamma * (tt)))^(k)
            }
            else if (type == "PGFInv.Der")
            {
              function(tt, gamma)
                -log(1 - gamma) / gamma * (1 - gamma)^(tt)
            }
          }
     }
     t@simul <- function(z, gamma) copula::rlog(z, gamma)
     t@theta <- vector("numeric")
     t@cop <- function(gamma, dim) Frank(gamma, dim)

     t
})
