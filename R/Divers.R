#' Automatic assigner
#'
#' Assign vector of values to variables automatically.
#' @param type Type of attribution: 0 = vector, 1 = row and 2 = column
#' @param data Data used (has to be the same as input)
#' @param output Desired variable for attribution
#' @param input Data variable: vector (0) or matrix (1 or 2)
#' @return Vector of attribution (expression)
#' @keywords internal
#' @export

var.attrib <- compiler::cmpfun(function(type, data, output = "x", input = "k", use = "param")
{
  output.length = length(output)

  if (type == 1) ## Row attrib.
  {
    dim <- length(data[,1])
    in1 <- paste(input, "[", 1:dim, ",]", sep = "")
    if (use != "param")
      out1 <- paste(output, 1:dim, sep = "")
    else
      out1 <- output[1:output.length]
  }

  if (type == 2) ## Column attrib.
  {
    dim <- length(data[1,])
    in1 <- paste(input, "[,", 1:dim, "]", sep = "")
    if (use != "param")
      out1 <- paste(output, 1:dim, sep = "")
    else
      out1 <- output[1:output.length]
  }

  if (type == 0)
  {
    dim <- length(data)
    in1 <- paste(input, "[", 1:dim, "]", sep = "")
    if (use != "param")
      out1 <- paste(output, 1:dim, sep = "")
    else
      out1 <- output[1:output.length]
  }

  expr <- "zzzzz <- zzzzzz"
  expr2 <- numeric(dim)
  for (i in 1:dim)
  {
    rep1 <- stringr::str_replace(expr, "zzzzz", out1[i])
    expr2[i] <- stringr::str_replace(rep1, "zzzzzz", in1[i])
  }
  parse(text = expr2)
})

#' MLE for copulas
#'
#' Estimate the parameters of a copula.
#' @param start Initial parameter
#' @param data Data used for the estimation
#' @param gradient (Not yet included)
#' @param low Lower bound for the research
#' @param up Upper bound fo the research
#' @param FUN The density of the copula to estimate (an expression)
#' @param param.type The parameters to estimate (character)
#' @importFrom stats nlminb
#' @return The estimation.
#' @seealso \code{\link{var.attrib}}
#' @export

maxvCopula <- compiler::cmpfun(function(start, data, low, up, FUN, gradient = NULL, param.type = "alpha")
{
  k <- data
  eval(var.attrib(2, k, input = "k", output = "x", use = "var"))
  xi <- paste("x", 1:length(k[1,]), sep = "")
  exp1 <- eval(parse(text = paste(xi, collapse = " * ")))


  logv <- function(param)
  {
    eval(var.attrib(0, param, input = "param", output = param.type))
    exp2 <- eval(FUN)
    -sum(log(exp1 * exp2))
  }

  nlminb(start, logv, lower = low, upper = up)
})

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

#' Extension of pairs
#'
#' Matrix of scatterplots with respective Kendall's tau
#' @param x Matrix of data
#' @param labels Vector of character (names of the variables)
#' @param cex Scale of the graph
#' @param cex.labels Numeric. The size of the labels' font
#' @param digits Number of digits to be shown for Kendall's tau
#' @return The pair graphs with respective Kendall tau
#' @importFrom graphics pairs par text
#' @importFrom Kendall Kendall
#' @seealso \code{\link{pairs}}
#' @export

pairs2 <- compiler::cmpfun(function(x, cex = 1, labels = paste("u", 1:length(x[1,]), sep = ""), cex.labels = 1, digits = 2)
{
  panel.cor <- function(x, y, cex.cor, ...)
  {
    usr <- par("usr"); on.exit(par(usr))
    par(usr = c(0, 1, 0, 1))
    # correlation coefficient
    r <- Kendall::Kendall(x, y)
    txt <- format(c(r, 0.123456789), digits = digits)[1]
    expr <- substitute(hat(tau) ==  u, list(u = txt))
    text(0.5, 0.6, expr, cex = 2)
  }
  pairs(x, upper.panel = panel.cor, labels = labels, cex.labels = cex.labels, cex = cex)
})

#' Shuffle of Min for the Frechet Lower Bound Uniforms
#'
#' Shuffle of min uniforms
#' @param n Number of simulations (numeric)
#' @param a Lower bounds (vector)
#' @param b Upper bounds (vector)
#' @param seed set.seed value
#' @param unif If true, returns the shuffles uniforms as well
#' @importFrom graphics lines
#' @importFrom graphics plot
#' @importFrom graphics points
#' @export


shuff <- compiler::cmpfun(function(n, a, b, seed = sample(1:200000, 1), unif = FALSE){

  set.seed(seed)
  par(pch = 16, cex = 0.8)
  x <- runif(n)
  y <- 1 - x

  x2.1 <- x[which({x >= a[1] & x <= b[1]})]

  x2.2 <- x[-which({x >= a[1] & x <= b[1]})]
  y2.2 <- 1 - x2.2
  res.x2 <- x2.1
  res.y2 <- 1 - res.x2 + runif(1, -min(1 - res.x2), 1 - max(1 - res.x2))

  res.x.ini <- c(res.x2)
  res.y.ini <- c(res.y2)

  if (unif == FALSE){
  if (length(a) == 1)
  {
    plot(x2.2, y2.2, type = "p")
    points(res.x2, res.y2, type = "p", xlim = c(0,1), ylim = c(0, 1), col = 2)
  }

  else
    plot(res.x2, res.y2, type = "p", xlim = c(0,1), ylim = c(0, 1), col = 2)

  #lines(x = c(min(res.x2), min(res.x2)), y = c(0, 1))
  lines(x = c(min(res.x2), min(res.x2)), y = c(min(1 - a[1], max(res.y2)), max(1 - a[1], max(res.y2))), lty = 3)
  lines(x = c(max(res.x2), max(res.x2)), y = c(min(1 - b[1], min(res.y2)), max(1 - b[1], min(res.y2))), lty = 3)}

  #lines(x = c(max(res.x2), max(res.x2)), y = c(0, 1))

  if (length(a) >= 2){

    for (i in 2:length(a))
    {
      x <- x2.2
      y <- 1 - x2.2

      x2.1 <- x[which({x >= a[i] & x <= b[i]})]

      x2.2 <- x[-which({x >= a[i] & x <= b[i]})]

      y2.2 <- 1 - x2.2
      res.x2 <- x2.1
      res.y2 <- 1 - res.x2 + runif(1, -min(1 - res.x2), 1 - max(1 - res.x2))

      res.x.ini <- c(res.x.ini, res.x2)
      res.y.ini <- c(res.y.ini, res.y2)

      if (unif == FALSE){
      if (i == length(a))
      {
        points(x2.2, y2.2, type = "p")
        points(res.x2, res.y2, type = "p", xlim = c(0,1), ylim = c(0, 1), col = 2)
      }
      else
        points(res.x2, res.y2, type = "p", xlim = c(0,1), ylim = c(0, 1), col = 2)

      #lines(x = c(min(res.x2), min(res.x2)), y = c(0, 1), lty = 3)
      lines(x = c(min(res.x2), min(res.x2)), y = c(min(1 - a[i], max(res.y2)), max(1 - a[i], max(res.y2))), lty = 3)
      lines(x = c(max(res.x2), max(res.x2)), y = c(min(1 - b[i], min(res.y2)), max(1 - b[i], min(res.y2))), lty = 3)}

      #lines(x = c(max(res.x2), max(res.x2)), y = c(0, 1), lty = 3)
    }}

  if (unif)
  {
    res.x.ini <- c(res.x.ini, x2.2)
    res.y.ini <- c(res.y.ini, y2.2)
    return(cbind(res.x.ini, res.y.ini))
  }
})

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
              LTheta = "character"),
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
              PM = "character"),
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
              PM = "character"),
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
              LTheta = "character"),
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
              LTheta = "character"),
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

