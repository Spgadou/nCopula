#' Copula Phi templates
#'
#' Copula Phi templates
#' @keywords internal
#' @export

gene.cop <- list("Clayton" = expression("phi" <- "exp(log((z) + 1)*(-1/alpha))",
                                        "phi.inv" <- "((z)^(-alpha) - 1)"),
                 "4.2.2" = expression("phi" <- "1 - (z)^(1/alpha)",
                                      "phi.inv" <- "(1 - (z))^alpha"),
                 "AMH" = expression("phi" <- "(1 - alpha) / (exp(z) - alpha)",
                                       "phi.inv" <- "log((1 - alpha*(1 - (z)))/(z))"),
                 "Gumbel" = expression("phi" <- "exp(-(z)^(1/alpha))",
                                       "phi.inv" <- "(-log(z))^alpha"),
                 "Frank" = expression("phi" <- "-log(exp(-(z)) * (exp(-alpha) - 1) + 1)/alpha",
                                      "phi.inv" <- "-log((exp(-alpha*(z)) - 1) / (exp(-alpha) - 1))"),
                 "4.2.6" = expression("phi" <- "1 - (1 + (z))^(1/alpha)",
                                      "phi.inv" <- "-log(1 - (1 - (z))^alpha)"),
                 "4.2.7" = expression("phi" <- "(exp(-(z)) + alpha - 1) / alpha",
                                      "phi.inv" <- "-log(alpha * (z) + (1 - alpha))"),
                 "4.2.13" = expression("phi" <- "exp(1 - (1 - (z))^(1/alpha))",
                                       "phi.inv" <- "(1 - log(z))^alpha"),
                 "4.2.15" = expression("phi" <- "((z)^(1/alpha) - 1)^alpha",
                                       "phi.inv" <- "(1 + (z)^(1/alpha))^alpha"),
                 "4.2.17" = expression("phi" <- "(exp(-(z)) * (2^(-alpha) - 1) + 1)^(-1/alpha) - 1",
                                       "phi.inv" <- "-log(((1 + (z))^(-alpha) - 1) / (2^(-alpha) - 1))"),
                 "4.2.20" = expression("phi" <- "(log(z + exp(1)))^(-1/alpha)",
                                       "phi.inv" <- "exp((t)^(-alpha)) - exp(1)"),
                 "4.2.21" = expression("phi" <- "1 - (1 - (1 - (z))^alpha)^(1/alpha)",
                                       "phi.inv" <- "1 - (1 - (1 - (z))^alpha)^(1/alpha)"))

#' Archimedean copulas family
#'
#' Archimedean copulas
#' @keywords internal
#' @exportClass archm

setClass("archm", list(phi = "character", phi.inv = "character", theta = "function", depend = "character", dimension = "numeric", parameter = "numeric", name = "character"),
         sealed = TRUE)

#' Clayton copula class
#'
#' Clayton copula
#' @keywords internal
#' @exportClass clayton


setClass("clayton",
         list(theta = "function",
              depend = "character",
              phi = "character",
              dens = "character",
              phi.inv = "character",
              rBiv = "function",
              dimension = "numeric", parameter = "numeric", name = "character"),
         contains = "archm", sealed = TRUE)


#' Frank copula class
#'
#' Frank copula
#' @keywords internal
#' @exportClass frank

setClass("frank",
         list(theta = "function", depend = "character",rBiv = "function",phi = "character", phi.inv = "character", dimension = "numeric", parameter = "numeric", name = "character"),
         contains = "archm", sealed = TRUE)

#' AMH copula class
#'
#' AMH copula
#' @keywords internal
#' @exportClass amh

setClass("amh",
         list(theta = "function", depend = "character",rBiv = "function",phi = "character", phi.inv = "character", dimension = "numeric", parameter = "numeric", name = "character"),
         contains = "archm", sealed = TRUE)

#' Gumvel copula class
#'
#' Gumbel copula
#' @keywords internal
#' @exportClass gumbel

setClass("gumbel",
         list(theta = "function", depend = "character",rBiv = "function",phi = "character", phi.inv = "character", dimension = "numeric", parameter = "numeric", name = "character"),
         contains = "archm", sealed = TRUE)


#' Mother Class
#'
#' @keywords internal
#' CompCop structure
#' @exportClass Mother

setClass("Mother", list(parameter = "numeric", structure = "list", arg = "numeric", dimension = "numeric"), sealed = TRUE)

#' Child Class
#'
#' @keywords internal
#' CompCop structure
#' @exportClass Child

setClass("Child", list(parameter = "numeric", arg = "numeric", dimension = "numeric"), sealed = TRUE)

#' Log-Child Class
#'
#' @keywords internal
#' CompCop structure
#' @exportClass Log_Child

setClass("Log_Child",
         list(name = "character",
              type = "character",
              dimension = "numeric",
              parameter = "numeric",
              arg = "numeric",
              obj = "character",
              Param = "character",
              Laplace = "character",
              LaplaceInv = "character",
              PGF = "character",
              PGFInv = "character",
              simul = "function",
              theta = "numeric",
              LTheta = "character",
              cop = "function"),
         contains = "Child", sealed = TRUE)

#' Log-Mother Class
#'
#' @keywords internal
#' CompCop structure
#' @exportClass Log_Mother

setClass("Log_Mother",
         list(name = "character",
              type = "character",
              dimension = "numeric",
              parameter = "numeric",
              arg = "numeric",
              structure = "list",
              obj = "character",
              Param = "character",
              Laplace = "character",
              LaplaceInv = "character",
              PGF = "character",
              PGFInv = "character",
              simul = "function",
              theta = "numeric",
              PM = "character",
              cop = "function"),
         contains = "Mother", sealed = TRUE)

#' Sibuya-Child Class
#'
#' @keywords internal
#' CompCop structure
#' @exportClass Sibuya_Child

setClass("Sibuya_Child",
         list(name = "character",
              type = "character",
              dimension = "numeric",
              parameter = "numeric",
              arg = "numeric",
              obj = "character",
              Param = "character",
              Laplace = "character",
              LaplaceInv = "character",
              PGF = "character",
              PGFInv = "character",
              simul = "function",
              theta = "numeric",
              LTheta = "character"),
         contains = "Child", sealed = TRUE)

#' Sibuya-Mother Class
#'
#' @keywords internal
#' CompCop structure
#' @exportClass Sibuya_Mother

setClass("Sibuya_Mother",
         list(name = "character",
              type = "character",
              dimension = "numeric",
              parameter = "numeric",
              arg = "numeric",
              structure = "list",
              obj = "character",
              Param = "character",
              Laplace = "character",
              LaplaceInv = "character",
              PGF = "character",
              PGFInv = "character",
              simul = "function",
              theta = "numeric",
              PM = "character"),
         contains = "Mother", sealed = TRUE)

#' Geo-Mother Class
#'
#' CompCop structure
#' @keywords internal
#' @exportClass Geo_Mother

setClass("Geo_Mother",
         list(name = "character",
              type = "character",
              dimension = "numeric",
              parameter = "numeric",
              structure = "list",
              arg = "numeric",
              obj = "character",
              Param = "character",
              Laplace = "character",
              LaplaceInv = "character",
              PGF = "character",
              PGFInv = "character",
              simul = "function",
              theta = "numeric",
              plot = "function",
              PM = "character",
              cop = "function",
              Der = "function",
              FUN = "function"),
         contains = "Mother", sealed = TRUE)

#' Geo-Child Class
#'
#' CompCop structure
#' @keywords internal
#' @exportClass Geo_Child

setClass("Geo_Child",
         list(name = "character",
              type = "character",
              dimension = "numeric",
              parameter = "numeric",
              arg = "numeric",
              obj = "character",
              Param = "character",
              Laplace = "character",
              LaplaceInv = "character",
              PGF = "character",
              PGFInv = "character",
              simul = "function",
              theta = "numeric",
              LTheta = "character",
              cop = "function",
              Der = "function",
              FUN = "function"),
         contains = "Child", sealed = TRUE)


#' Gamma-Child Class
#'
#' CompCop structure
#' @keywords internal
#' @exportClass Gamma_Child

setClass("Gamma_Child",
         list(name = "character",
              type = "character",
              dimension = "numeric",
              parameter = "numeric",
              arg = "numeric",
              obj = "character",
              Param = "character",
              Laplace = "character",
              LaplaceInv = "character",
              PGF = "character",
              PGFInv = "character",
              simul = "function",
              theta = "numeric",
              LTheta = "character",
              Der = "function",
              FUN = "function"),
         contains = "Child", sealed = TRUE)


#' Show Method for Copulas
#'
#' @keywords internal
#' Show

setMethod("show",
          "archm",
          definition = function(object)
          {
            cat(object@name, "\n")
            cat("\n")
            cat("   Dimension :", object@dimension, "\n")
            cat("   Parameter :", object@parameter, "\n")
            cat("\n")
          }
)

#' Show Method Function for Compounding
#'
#' @param object S4 object
#' @param indent Initial spacing for children
#' @param delta.ident Added space for each iteration
#' @param label Internal use (must not be changed)
#' @keywords internal
#' @export

Moth <- compiler::cmpfun(function(object, indent = "", delta.indent = 3, label = NA){

  mkBlanks <- function(n) paste(rep.int(" ", n), collapse = "")
  space <- mkBlanks(nIS <- nchar(indent))

  if (object@type == "Mother")
  {

    nk <- paste0(seq_len(length(object@structure)), ")")
    if (is.na(label))
      cat(space, object@name, "\n")
    else
      cat(space, label, object@name, "\n")
    cat("\n")
    cat("  ", space, "Type     :", object@type, "\n")
    cat("  ", space, "Parameter:", object@parameter, "\n")
    cat("  ", space, "Dimension:", object@dimension, "\n")
    if (is.null(object@arg) == FALSE)
      cat("  ", space, "Arguments:", object@arg, "\n")
    cat("  ", space, "Children:", "\n")
    cat("\n")

    space <- mkBlanks(nIS + delta.indent)
    for (k in 1:length(object@structure))
    {
      Moth(object@structure[[k]], indent = paste0(space), label = nk[k])
    }
  }

  else
  {
    cat(space, label, object@name, "\n")
    cat("\n")
    cat("  ", space, "Type     :", object@type, "\n")
    cat("  ", space, "Parameter:", object@parameter, "\n")
    cat("  ", space, "Dimension:", object@dimension, "\n")
    cat("  ", space, "Arguments:", object@arg, "\n")
    cat("\n")
  }
})

#' @rdname Moth
#' @export

Moth2 <- compiler::cmpfun(function(object, indent = "", delta.indent = 3, label = NA){

  mkBlanks <- function(n) paste(rep.int(" ", n), collapse = "")
  space <- mkBlanks(nIS <- nchar(indent))

  if (object@type == "Mother")
  {
    chil <- object@dimension - (length(object@arg) > 1 || {length(object@arg) == 1 && object@arg != 0}) * length(object@arg)

    if ((length(object@arg) > 1 || {length(object@arg) == 1 && object@arg != 0}))
      ui <- paste("u", object@arg, sep = "")

    nk <- paste0(seq_len(length(object@structure)), ")")
    if (is.na(label))
      cat(space, paste0(object@name, ":"), paste0(object@dimension + length(object@arg) * (length(object@arg) > 1 || {length(object@arg) == 1 && object@arg != 0}), "-dimensional "))
    else
      cat(space,label,paste0(object@name, ":"), paste0(object@dimension, "-dimensional "))

    cat(paste0("'", object@type, "'"), "function", "with parameter", round(object@parameter, 4))
    if (object@type == "Mother")
    {
      cat("\n")
      cat(space)
      cat(space, "composed of", paste0("(", chil, ")"), paste0("child", if (chil > 1){ "ren"}))

      if ((length(object@arg) > 1 || {length(object@arg) == 1 && object@arg != 0}))
        cat(" and", paste0(paste0("(",paste(ui, collapse = ", "), ")"), ":"), "\n")
      else
        cat(":", "\n")
    }
    else
      cat(":", "\n")
    cat("\n")

    space <- mkBlanks(nIS + delta.indent)
    for (k in 1:length(object@structure))
    {
      Moth2(object@structure[[k]], indent = paste0(space), label = nk[k])
    }
  }

  else
  {
    ui <- paste("u", object@arg, sep = "")

    cat(space, label, paste0(object@name, ":"), "")
    cat(paste0(object@dimension, "-dimensional "))
    cat(paste0("'", object@type, "'"), "function", "with parameter", round(object@parameter, 4), "\n")
    cat(space, "  ", "composed of", paste0("(",paste(ui, collapse = ", "), ")"), "\n")
    cat("\n")
  }
})

#' @rdname Moth
#' @export

Chil <- compiler::cmpfun(function(object, space = ""){

    ui <- paste("u", object@arg, sep = "")

    cat(space, paste0(object@name, ":"), "")
    cat(paste0(object@dimension, "-dimensional "))
    cat(paste0("'", object@type, "'"), "function", "with parameter", round(object@parameter, 4), "\n")
    cat(space, "composed of", paste0("(",paste(ui, collapse = ", "), ")"), "\n")
  }
)

#' Show Method for Compounding Functions - Mother Class
#'
#' @keywords internal
#' Show

setMethod("show",
          "Mother",
          definition = function(object) Moth2(object))

#' Show Method for Compounding Functions - Child Class
#'
#' @keywords internal
#' Show

setMethod("show",
          "Child",
          definition = function(object) Chil(object))

