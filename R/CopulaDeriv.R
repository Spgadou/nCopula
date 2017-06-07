#' Construction of a Mother or Child Class Object
#'
#' @description Constructs either a Mother or Child class object for
#' a given parameter, arguments, and nesting structure.
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

#' @rdname LOG
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

#' @rdname LOG
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

#' @rdname LOG
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

#' Construction of an Archimedean Copula Class Object
#'
#' @description Constructs an archm class object with a given parameter and
#' dimension.
#' @param dim Dimension of the copula (>= 2), which is, by default, 2
#' @param param Parameter of the copula
#' @importFrom copula rlog
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
  dep.param <- "1/alpha"
  rBiv <- function(n, alpha, u) (u^(-alpha) * (runif(n)^(-alpha / (alpha + 1)) - 1) + 1)^(-1/alpha)
  th <- function(z, alpha) rgamma(z, alpha, 1)

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

    expr2 <- stringr::str_replace_all(tt@Laplace, "z", nu)
    densit <- paste("(", expr1, ") * (", expr2, ")", sep = "")
    densit <- stringr::str_replace_all(densit, "alpha", "(1/alpha)")

    new("clayton",
        phi = phi,
        phi.inv = phi.inv,
        rBiv = rBiv,
        theta = th,
        depend = dep.param,
        dimension = dim,
        parameter = param,
        dens = densit,
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
        name = "Clayton copula")
  }
})

#' @rdname Clayton
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

#' @rdname Clayton
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

#' @rdname Clayton
#' @export

Gumbel <- compiler::cmpfun(function(param, dim = 2L)
{
  if (param < 1)
    stop("Wrong 'param' input")

  verif <- eval(parse(text = paste0(dim, "L")))

  if (!is.integer(verif) || dim <= 1)
    stop("The dimension must be an integer greater than or equal to 2")

  phi <- "exp(-(z)^(1/alpha))"
  phi.inv <- "(-log(z))^alpha"
  dep.param <- "alpha"
  th <- function(z, alpha) copula::rstable1(z, 1/alpha, 1, cos(pi/(2*alpha))^alpha, 0, 1)

  new("gumbel",
      phi = phi,
      phi.inv = phi.inv,
      theta = th,
      depend = dep.param,
      dimension = dim,
      parameter = param,
      name = "Gumbel copula")
})

#' Density, Cdf, and Random Number Generator for Copulas
#'
#' @param n Number of realisations
#' @param level Imbrication levels
#' @param max.dim Safe upper bound for final dimension (2000 by default)
#' @param copula An Archimedean copula class object
#' @param vector If false, returns a function with (x_1, x_2, ..., x_dim, alpha) as arguments.
#' @param express If true, returns an expression.
#' @param code If true, copies the LaTeX code to clipboard.
#' @param operator Type of cumputer used (only necessary in the case of internal problem)
#' @return Either an expression, function, code, or sampled data.
#' @details rCop2 combines both the features of rCop and rCompCop in one, allowing the input
#' of Archimedean copula classes, Mother classes, and Child classes.
#' @author Simon-Pierre Gadoury
#' @export

dCop <- compiler::cmpfun(function(copula, vector = TRUE, express = FALSE, code = FALSE, operator = "")
{

  system <- Sys.info()["sysname"]

  phi <- copula@phi
  dim <- copula@dimension
  phi.inv <- copula@phi.inv

  ui <- paste("x", 1:dim, sep = "")
  ui2 <- paste("x[", 1:dim, "]", sep = "")
  expr1 <- stringr::str_replace(phi.inv, "z", ui)
  expr1 <- paste(expr1, collapse = " + ")
  expr2 <- stringr::str_replace(phi, "z", expr1)

  der <- list()
  der[[1]] <- Deriv::Deriv(expr2, ui[1], cache.exp = FALSE)
  for (i in 2:dim)
    der[[i]] <- Deriv::Deriv(der[[i - 1]], ui[i], cache.exp = FALSE)
  der <- der[[dim]]

  if (vector == FALSE && express == FALSE)
  {
    res <- stringr::str_replace_all(der, "dim", dim)

    t1 <- "function(z)"
    t2 <- paste(ui, collapse = ", ")
    t3 <- paste(c(t2, "alpha"), collapse = ", ")
    input <- stringr::str_replace(t1, "z", t3)
    input2 <- paste(c(input, res), collapse = " ")
    res2 <- parse(text = input2)
    return(eval(res2))
  }

  if (vector == FALSE && express == TRUE)
  {
    res <- stringr::str_replace_all(der, "dim", dim)
    return(parse(text = res))
  }

  for (i in 1:dim)
  {
    j <- dim - i + 1
    der <- stringr::str_replace_all(der, ui[j], ui2[j])
  }

  res <- stringr::str_replace_all(der, "dim", dim)

  if (express)
    return(parse(text = res))

  if (code)
  {
    num2 <- paste("_", 1:dim, sep = "")
    num2 <- paste(c(rep("u", dim), num2), collapse = "")
    for (i in 1:dim)
    {
      j <- dim - i + 1
      res <- stringr::str_replace_all(res, ui2[j], num2[j])
    }

    tex <- paste("$", Ryacas::yacas(Ryacas::TeXForm(res), retclass = "unquote"), "$", sep = "")

    if (system == "Darwin" || operator == "Mac"){
      clip <- pipe("pbcopy", "w")
      write(tex, file = clip)
      close(clip)
      return("Code succesfully copied to clipboard")}

    if (system == "Windows" || operator == "Windows" || operator == "Linux" || system == "Linux")
    {
      writeClipboard(tex)
      return("Code succesfully copied to clipboard")
    }

    else
      stop("Due to internal error, manually use the 'operator' argument: Windows -> 'Windows', OSX -> 'Darwin, Linux -> 'Linux'")
  }

  else
  {
    res2 <- function(x, alpha) eval(parse(text = res))
    return(res2)
  }
})

#' @rdname dCop
#' @export

pCop <- compiler::cmpfun(function(copula, vector = TRUE, express = FALSE, code = FALSE, operator = "")
{

  system <- Sys.info()["sysname"]

  phi <- copula@phi
  dim <- copula@dimension
  phi.inv <- copula@phi.inv

  ui <- paste("x", 1:dim, sep = "")
  ui2 <- paste("x[", 1:dim, "]", sep = "")
  expr1 <- stringr::str_replace_all(phi.inv, "z", ui)
  expr1 <- paste(expr1, collapse = " + ")
  expr2 <- stringr::str_replace_all(phi, "z", expr1)

  if (vector == FALSE && express == TRUE)
  {
    expr2 <- Deriv::Simplify(expr2)
    res <- stringr::str_replace_all(expr2, "dim", dim)
    return(parse(text = res))
  }

  if (vector == FALSE && express == FALSE)
  {
    expr2 <- Deriv::Simplify(expr2)
    res <- stringr::str_replace_all(expr2, "dim", dim)

    t1 <- "function(z)"
    t2 <- paste(ui, collapse = ", ")
    t3 <- paste(c(t2, "alpha"), collapse = ", ")
    input <- stringr::str_replace_all(t1, "z", t3)
    input2 <- paste(c(input, res), collapse = " ")
    res2 <- parse(text = input2)
    return(eval(res2))
  }

  for (i in 1:dim)
  {
    j <- dim - i + 1
    expr2 <- stringr::str_replace_all(expr2, ui[j], ui2[j])
  }

  expr2 <- stringr::str_replace_all(expr2, "dim", dim)
  expr2 <- Deriv::Simplify(expr2)

  if (express)
  {
    return(parse(text = expr2))
  }

  if (code)
  {
    res <- expr2

    num2 <- paste("_", 1:dim, sep = "")
    num2 <- paste(c(rep("u", dim), num2), collapse = "")
    for (i in 1:dim)
    {
      j <- dim - i + 1
      expr2 <- stringr::str_replace_all(expr2, ui2[j], num2[j])
    }

    tex <- paste("$", Ryacas::yacas(Ryacas::TeXForm(expr2), retclass = "unquote"), "$", sep = "")

    if (system == "Darwin" || operator == "Mac"){
      clip <- pipe("pbcopy", "w")
      write(tex, file = clip)
      close(clip)
      return("Code succesfully copied to clipboard")}

    if (system == "Windows" || operator == "Windows" || operator == "Linux" || system == "Linux")
    {
      writeClipboard(tex)
      return("Code succesfully copied to clipboard")
    }

    else
      stop("Due to internal error, manually change the 'operator' argument: Windows -> 'Windows', OSX -> 'Darwin, Linux -> 'Linux'")
  }

  else
  {
    t1 <- "function(z)"
    t3 <- paste(c("x", "alpha"), collapse = ", ")
    input <- stringr::str_replace_all(t1, "z", t3)
    input2 <- paste(c(input, expr2), collapse = " ")
    res2 <- parse(text = input2)
    return(eval(res2))
  }
})

#' @rdname dCop
#' @export

rCop <- compiler::cmpfun(function(n, copula)
{
  sim <- copula@theta
  param <- copula@parameter
  dim <- copula@dimension

  param2 <- copula@depend
  param2 <- stringr::str_replace_all(param2, "alpha", param)
  param2 <- eval(parse(text = param2))

  phi <- stringr::str_replace_all(copula@phi, "alpha", param)
  phi <- parse(text = stringr::str_replace_all(phi, "z", "res"))

  res <- t(-log(runif(dim*n)) / matrix(sim(n, param2), ncol = n, nrow = dim, byrow = TRUE))

  eval(phi)
})

#' @rdname dCop
#' @export

rCop2 <- compiler::cmpfun(function(n, copula, level = "", max.dim = 2000)
{
  tt <- sum(c(names(getClass("Mother")@subclasses),names(getClass("Child")@subclasses))  == class(copula)[1])

  if (tt == 0){

    if (copula@dimension == 2)
    {
      u1 <- runif(n)
      u2 <- copula@rBiv(n, copula@parameter, u1)
      cbind(u1, u2)
    }

    else
    {
    sim <- copula@theta
    param <- copula@parameter
    dim <- copula@dimension

    param2 <- copula@depend
    param2 <- stringr::str_replace_all(param2, "alpha", param)
    param2 <- eval(parse(text = param2))

    phi <- stringr::str_replace_all(copula@phi, "alpha", param)
    phi <- parse(text = stringr::str_replace_all(phi, "z", "res"))

    res <- t(-log(runif(dim*n)) / matrix(sim(n, param2), ncol = n, nrow = dim, byrow = TRUE))

    eval(phi)}
  }

  else
  {
    FUN <- copula

    if (FUN@type == "Child")
    {
      FUN@theta <- FUN@simul(n, FUN@parameter)
      laplace <- stringr::str_replace_all(FUN@Laplace, FUN@Param, FUN@parameter)
      th <- "-log(runif(FUN@dimension * n)) / FUN@theta"
      laplace <- parse(text = stringr::str_replace_all(laplace, "z", th))
      return(matrix(eval(laplace), ncol = FUN@dimension, nrow = n))
    }

    t <- FUN
    res <- vector("list", max.dim)
    FUN <- list()

    th.pos1 <- list()
    m.pos1 <- list()

    th.pos1[[1]] <- list()
    m.pos1[[1]] <- list()

    FUN[[1]] <- list()
    FUN[[1]][[1]] <- t


    for (i in 1:level)
    {
      FUN[[i + 1]] <- list()

      if (i == 1){

        FUN[[1]][[1]]@theta <- FUN[[1]][[1]]@simul(n, FUN[[1]][[1]]@parameter)

        FUN[[1]][[1]]@PGF <- stringr::str_replace_all(FUN[[1]][[1]]@PGF, FUN[[1]][[1]]@Param, FUN[[1]][[1]]@parameter)

        if (FUN[[1]][[1]]@type == "Mother"){
          typ <- numeric(FUN[[1]][[1]]@dimension)
          for (j in 1:length(FUN[[1]][[1]]@structure))
            typ[j] <- FUN[[1]][[1]]@structure[[j]]@type}

        th.pos1[[1]][[1]] <- which(typ == "Child")
        m.pos1[[1]][[1]] <- which(typ == "Mother")

        if (length(FUN[[1]][[1]]@arg) > 1 || {length(FUN[[1]][[1]]@arg) == 1 && FUN[[1]][[1]]@arg != 0})
        {
          for (j in 1:length(FUN[[1]][[1]]@arg))
          {
            laplace <- stringr::str_replace_all(FUN[[1]][[1]]@Laplace, FUN[[1]][[1]]@Param, FUN[[1]][[1]]@parameter)
            th <- "-log(runif(n)) / FUN[[1]][[1]]@theta"
            laplace <- parse(text = stringr::str_replace_all(laplace, "z", th))

            if (!is.null(res[[FUN[[1]][[1]]@arg[j]]]))
              warning("Multiply defined 'unif' arguments will result in the last one processed")

            res[[FUN[[1]][[1]]@arg[j]]] <- eval(laplace)
          }
        }

        if (length(th.pos1[[1]][[1]]) != 0){
          for (j in 1:length(th.pos1[[1]][[1]]))
          {
            laplace <- stringr::str_replace_all(FUN[[1]][[1]]@structure[[th.pos1[[1]][[1]][j]]]@Laplace, FUN[[1]][[1]]@structure[[th.pos1[[1]][[1]][j]]]@Param, FUN[[1]][[1]]@structure[[th.pos1[[1]][[1]][j]]]@parameter)
            fbarre <- stringr::str_replace_all(FUN[[1]][[1]]@PGF, "z", laplace)

            argg <- FUN[[1]][[1]]@structure[[th.pos1[[1]][[1]][j]]]@arg

            th2 <- vapply(FUN[[1]][[1]]@theta, function(x) sum(FUN[[1]][[1]]@structure[[th.pos1[[1]][[1]][j]]]@simul(x, FUN[[1]][[1]]@structure[[th.pos1[[1]][[1]][j]]]@parameter)), rep(0, 1))

            res2 <- list()
            for (l in 1:length(argg))
            {
              th <- -log(runif(n)) / th2

              fbarre <- stringr::str_replace_all(fbarre, "z", "th")

              if (!is.null(res[[argg[l]]]))
                warning("Multiply defined 'unif' arguments will result in the last one processed")

              res[[argg[l]]] <- eval(parse(text = fbarre))
            }
          }}

        if (length(m.pos1[[1]][[1]]) != 0){
          for (j in 1:length(m.pos1[[1]][[1]]))
          {

            FUN[[1]][[1]]@structure[[m.pos1[[1]][[1]][j]]]@PGF <- stringr::str_replace_all(FUN[[1]][[1]]@structure[[m.pos1[[1]][[1]][j]]]@PGF, FUN[[1]][[1]]@structure[[m.pos1[[1]][[1]][j]]]@Param, FUN[[1]][[1]]@structure[[m.pos1[[1]][[1]][j]]]@parameter)

            FUN[[1]][[1]]@structure[[m.pos1[[1]][[1]][j]]]@PGF <- stringr::str_replace_all(FUN[[1]][[1]]@PGF, "z", FUN[[1]][[1]]@structure[[m.pos1[[1]][[1]][j]]]@PGF)

            FUN[[1]][[1]]@structure[[m.pos1[[1]][[1]][j]]]@theta <- vapply(FUN[[1]][[1]]@theta, function(x) sum(FUN[[1]][[1]]@structure[[m.pos1[[1]][[1]][j]]]@simul(x, FUN[[1]][[1]]@structure[[m.pos1[[1]][[1]][j]]]@parameter)), rep(0, 1))

            FUN[[i + 1]][[j]] <- FUN[[1]][[1]]@structure[[m.pos1[[1]][[1]][j]]]
          }}

      }

      else
      {
        th.pos1[[i]] <- list()
        m.pos1[[i]] <- list()

        for (k in 1:length(m.pos1[[i - 1]][[1]]))
        {

          typ <- numeric(FUN[[i]][[k]]@dimension)
          for (j in 1:length(FUN[[i]][[k]]@structure))
            typ[j] <- FUN[[i]][[k]]@structure[[j]]@type

          th.pos1[[i]][[k]] <- which(typ == "Child")
          m.pos1[[i]][[k]] <- which(typ == "Mother")

          if (FUN[[i]][[k]]@arg != 0)
          {
            for (j in 1:length(FUN[[i]][[k]]@arg))
            {
              laplace <- stringr::str_replace_all(FUN[[i]][[k]]@Laplace, FUN[[i]][[k]]@Param, FUN[[i]][[k]]@parameter)
              th <- "-log(runif(n)) / FUN[[i]][[k]]@theta"
              laplace <- parse(text = stringr::str_replace_all(laplace, "z", th))

              if (!is.null(res[[FUN[[i]][[k]]@arg[j]]]))
                warning("Multiply defined 'unif' arguments will result in the last one processed")

              res[[FUN[[i]][[k]]@arg[j]]] <- eval(laplace)
            }
          }

          if (length(th.pos1[[i]][[k]]) != 0){
            for (j in 1:length(th.pos1[[i]][[k]]))
            {
              laplace <- stringr::str_replace_all(FUN[[i]][[k]]@structure[[th.pos1[[i]][[k]][j]]]@Laplace, FUN[[i]][[k]]@structure[[th.pos1[[i]][[k]][j]]]@Param, FUN[[i]][[k]]@structure[[th.pos1[[i]][[k]][j]]]@parameter)

              fbarre <- stringr::str_replace_all(FUN[[i]][[k]]@PGF, "z", laplace)

              argg <- FUN[[i]][[k]]@structure[[th.pos1[[i]][[1]][j]]]@arg

              th2 <- vapply(FUN[[i]][[k]]@theta, function(x) sum(FUN[[i]][[k]]@structure[[th.pos1[[i]][[k]][j]]]@simul(x, FUN[[i]][[k]]@structure[[th.pos1[[i]][[k]][j]]]@parameter)), rep(0, 1))

              res2 <- list()
              for (l in 1:length(argg))
              {
                th <- -log(runif(n)) / th2

                fbarre <- parse(text = stringr::str_replace_all(fbarre, "z", "th"))

                if (!is.null(res[[argg[l]]]))
                  warning("Multiply defined 'unif' arguments will result in the last one processed")

                res[[argg[l]]] <- eval(fbarre)
              }

              #res[[argg[1]]] <- matrix(unlist(res2), ncol = FUN[[1]][[1]]@structure[[th.pos1[[1]][[1]][j]]]@dimension, nrow = n)
            }}

          if (length(m.pos1[[i]][[k]]) != 0){
            for (j in 1:length(m.pos1[[i]][[k]]))
            {

              FUN[[i]][[k]]@structure[[m.pos1[[i]][[k]][j]]]@PGF <- stringr::str_replace_all(FUN[[i]][[k]]@structure[[m.pos1[[i]][[k]][j]]]@PGF,
                                                                                             FUN[[i]][[k]]@structure[[m.pos1[[i]][[k]][j]]]@Param,
                                                                                             FUN[[i]][[k]]@structure[[m.pos1[[i]][[k]][j]]]@parameter)

              FUN[[i]][[k]]@structure[[m.pos1[[i]][[k]][j]]]@PGF <- stringr::str_replace_all(FUN[[i]][[k]]@PGF, "z", FUN[[i]][[k]]@structure[[m.pos1[[i]][[k]][j]]]@PGF)

              FUN[[i]][[k]]@structure[[m.pos1[[i]][[k]][j]]]@theta <- vapply(FUN[[i]][[k]]@theta, function(x) sum(FUN[[i]][[k]]@structure[[m.pos1[[i]][[k]][j]]]@simul(x, FUN[[i]][[k]]@structure[[m.pos1[[i]][[k]][j]]]@parameter)), rep(0, 1))

              FUN[[i + 1]][[j]] <- FUN[[i]][[k]]@structure[[m.pos1[[i]][[k]][j]]]

            }}
        }
      }

    }
    do.call(cbind, res)
  }
})



#' Density, Cdf, and Random Number Generator for Copulas Constructed Through Compounding
#'
#' @param n Number of realisations
#' @param level Number of imbrications
#' @param FUN Object of class Mother
#' @param func If true, returns a function
#' @param code If true, copies the LaTeX code to the clipboard
#' @param operator Type of cumputer used (only necessary in the case of internal problem)
#' @details rCop2 merges rCop and rCompCop2 in one function.
#' @return The copula
#' @export

dCompCop <- compiler::cmpfun(function(FUN, func = FALSE, code = FALSE, operator = ""){

  system <- Sys.info()["sysname"]

  dim.cop <- FUN@dimension
  PGFInv <- FUN@PGFInv

  num <- paste(seq(1, 2*dim.cop - 1, 2), seq(2, 2*dim.cop, 2), sep = "")

  pp <- list()
  pp2 <- list()
  pp3 <- list()
  org <- vector("numeric")

  u <- list()

  for (k in 1:dim.cop)
  {
    org[k] <- length(FUN@structure[[k]]@arg)

    ui <- paste("u", FUN@structure[[k]]@arg[1]:FUN@structure[[k]]@arg[org[k]], sep = "")

    u[[k]] <- ui

    param <- FUN@structure[[k]]@Param
    Laplace <- FUN@structure[[k]]@Laplace
    LaplaceInv <- FUN@structure[[k]]@LaplaceInv

    p <- numeric(length(param))
    for (j in 1:length(param))
    {
      p[j] <- paste(c(param[j], num[k]), collapse = "")
      Laplace <- stringr::str_replace_all(Laplace, param[j], p[j])
      LaplaceInv <- stringr::str_replace_all(LaplaceInv, param[j], p[j])
    }

    LaplaceInv<- stringr::str_replace_all(LaplaceInv, "z", PGFInv)

    inp <- list()
    for (i in 1:org[k])
      inp[[i]] <- stringr::str_replace_all(LaplaceInv, "z", ui[i])

    t <- paste(unlist(inp), collapse = " + ")

    pp2[[k]] <- stringr::str_replace_all(Laplace, "z", t)
  }

  input <- paste(unlist(pp2), collapse = " * ")

  PGF <- FUN@PGF
  res <- stringr::str_replace_all(PGF, "z", input)

  der <- list()
  ui <- unlist(u)
  der[[1]] <- Deriv::Deriv(res, "u1", cache.exp = FALSE)
  for (i in 2:length(ui))
    der[[i]] <- Deriv::Deriv(der[[i - 1]], ui[i], cache.exp = FALSE)
  res <- der[[length(ui)]]

  if (code)
  {
    num1 <- paste("_", seq(1, 2*dim.cop - 1, 2), seq(2, 2*dim.cop, 2), sep = "")
    for (i in 1:length(num))
      res <- stringr::str_replace_all(res, num[i], num1[i])
    tex <- paste("$", Ryacas::yacas(Ryacas::TeXForm(res), retclass = "unquote"), "$", sep = "")

    num2 <- paste("{", seq(1, 2*dim.cop - 1, 2), seq(2, 2*dim.cop, 2), "}", sep = "")
    for (i in 1:length(num))
      tex <- stringr::str_replace_all(tex, num[i], num2[i])

    if (system == "Darwin" || operator == "Mac"){
      clip <- pipe("pbcopy", "w")
      write(tex, file = clip)
      close(clip)
      return("Code succesfully copied to clipboard")}

    if (system == "Windows" || operator == "Windows" || operator == "Linux" || system == "Linux")
    {
      writeClipboard(tex)
      return("Code succesfully copied to clipboard")
    }

    else
      stop("Due to internal error, use the 'operator' argument: Windows -> 'Windows', OSX -> 'Darwin, Linux -> 'Linux'")
  }

  if (func)
  {
    t1 <- "function(z)"
    t2.1 <- paste(ui, collapse = ", ")
    t2.2 <- paste(unlist(pp3), collapse = ", ")
    t3 <- paste(c(t2.1, t2.2), collapse = ", ")
    input <- stringr::str_replace(t1, "z", t3)
    input2 <- paste(c(input, res), collapse = " ")
    res2 <- parse(text = input2)
    eval(res2)
  }

  else
  {
    return(res)
  }
})

#' @rdname dCompCop
#' @export

pCompCop <- compiler::cmpfun(function(FUN, func = FALSE, code = FALSE, operator = ""){

  system <- Sys.info()["sysname"]

  dim.cop <- FUN@dimension
  PGFInv <- FUN@PGFInv

  num <- paste(seq(1, 2*dim.cop - 1, 2), seq(2, 2*dim.cop, 2), sep = "")

  pp <- list()
  pp2 <- list()
  pp3 <- list()
  org <- vector("numeric")

  for (k in 1:dim.cop)
  {
    org[k] <- length(FUN@structure[[k]]@arg)

    ui <- paste("u", FUN@structure[[k]]@arg[1]:FUN@structure[[k]]@arg[org[k]], sep = "")

    param <- FUN@structure[[k]]@Param
    Laplace <- FUN@structure[[k]]@Laplace
    LaplaceInv <- FUN@structure[[k]]@LaplaceInv

    p <- numeric(length(param))
    for (j in 1:length(param))
    {
      p[j] <- paste(c(param[j], num[k]), collapse = "")
      Laplace <- stringr::str_replace_all(Laplace, param[j], p[j])
      LaplaceInv <- stringr::str_replace_all(LaplaceInv, param[j], p[j])
    }
    pp3[[k+1]] <- p

    LaplaceInv<- stringr::str_replace_all(LaplaceInv, "z", PGFInv)

    inp <- list()
    for (i in 1:org[k])
      inp[[i]] <- stringr::str_replace_all(LaplaceInv, "z", ui[i])

    t <- paste(unlist(inp), collapse = " + ")
    pp2[[k]] <- stringr::str_replace_all(Laplace, "z", t)
  }

  input <- paste(unlist(pp2), collapse = " * ")

  PGF <- FUN@PGF
  res <- stringr::str_replace_all(PGF, "z", input)
  #res <- Deriv::Simplify(res) ## Résultat simplifié ou non ?

  if (code)
  {
    num1 <- paste("_", seq(1, 2*dim.cop - 1, 2), seq(2, 2*dim.cop, 2), sep = "")
    for (i in 1:length(num))
      res <- stringr::str_replace_all(res, num[i], num1[i])
    tex <- paste("$", Ryacas::yacas(Ryacas::TeXForm(res), retclass = "unquote"), "$", sep = "")

    num2 <- paste("{", seq(1, 2*dim.cop - 1, 2), seq(2, 2*dim.cop, 2), "}", sep = "")
    for (i in 1:length(num))
      tex <- stringr::str_replace_all(tex, num[i], num2[i])

    if (system == "Darwin" || operator == "Mac"){
      clip <- pipe("pbcopy", "w")
      write(tex, file = clip)
      close(clip)
      return("Code succesfully copied to clipboard")}

    if (system == "Windows" || operator == "Windows" || operator == "Linux" || system == "Linux")
    {
      writeClipboard(tex)
      return("Code succesfully copied to clipboard")
    }

    else
      stop("Due to internal error, use the 'operator' argument: Windows -> 'Windows', OSX -> 'Darwin, Linux -> 'Linux'")
  }

  if (func)
  {
    t1 <- "function(z)"
    t2.1 <- paste(ui, collapse = ", ")
    t2.2 <- paste(unlist(pp3), collapse = ", ")
    t3 <- paste(c(t2.1, t2.2), collapse = ", ")
    input <- stringr::str_replace(t1, "z", t3)
    input2 <- paste(c(input, res), collapse = " ")
    res2 <- parse(text = input2)
    eval(res2)
  }

  else
  {
    return(res)
  }
})

#' @rdname dCompCop
#' @export

rCompCop <- compiler::cmpfun(function(n, FUN, level){

  t <- FUN
  res <- list()
  FUN <- list()

  th.pos1 <- list()
  m.pos1 <- list()

  th.pos1[[1]] <- list()
  m.pos1[[1]] <- list()

  FUN[[1]] <- list()
  FUN[[1]][[1]] <- t

  num1 <- list()

  for (i in 1:level)
  {
    FUN[[i + 1]] <- list()

    if (i == 1){

      FUN[[1]][[1]]@theta <- FUN[[1]][[1]]@simul(n, FUN[[1]][[1]]@parameter)

      typ <- numeric(FUN[[1]][[1]]@dimension)
      for (j in 1:length(FUN[[1]][[1]]@structure))
        typ[j] <- FUN[[1]][[1]]@structure[[j]]@type

      th.pos1[[1]][[1]] <- which(typ == "Child")
      m.pos1[[1]][[1]] <- which(typ == "Mother")

      num1[[1]] <- paste(1, 1, 1:FUN[[1]][[1]]@dimension, sep = "")

      if (length(FUN[[1]][[1]]@arg) > 1 || {length(FUN[[1]][[1]]@arg) == 1 && FUN[[1]][[1]]@arg != 0})
      {
        for (j in 1:length(FUN[[1]][[1]]@arg))
        {
          laplace <- stringr::str_replace_all(FUN[[1]][[1]]@Laplace, FUN[[1]][[1]]@Param, FUN[[1]][[1]]@parameter)
          th <- "-log(runif(n)) / FUN[[1]][[1]]@theta"
          laplace <- parse(text = stringr::str_replace_all(laplace, "z", th))
          res[[FUN[[1]][[1]]@arg[j]]] <- eval(laplace)
        }
      }

      for (j in 1:length(th.pos1[[1]][[1]]))
      {
        param <- paste(FUN[[1]][[1]]@structure[[th.pos1[[1]][[1]][j]]]@Param, num1[[1]][th.pos1[[1]][[1]][j]], sep = "")
        laplace <- stringr::str_replace_all(FUN[[1]][[1]]@structure[[th.pos1[[1]][[1]][j]]]@Laplace, FUN[[1]][[1]]@structure[[th.pos1[[1]][[1]][j]]]@Param, param)
        fbarre <- stringr::str_replace_all(FUN[[1]][[1]]@PGF, "z", laplace)

        fbarre <- stringr::str_replace_all(fbarre, param, FUN[[1]][[1]]@structure[[th.pos1[[1]][[1]][j]]]@parameter)
        fbarre <- stringr::str_replace_all(fbarre, FUN[[1]][[1]]@Param, FUN[[1]][[1]]@parameter)


        th <- matrix(-log(runif(FUN[[1]][[1]]@structure[[th.pos1[[1]][[1]][j]]]@dimension * n)) / vapply(FUN[[1]][[1]]@theta, function(x) sum(FUN[[1]][[1]]@structure[[th.pos1[[1]][[1]][j]]]@simul(x, FUN[[1]][[1]]@structure[[th.pos1[[1]][[1]][j]]]@parameter)), rep(0, 1)), ncol = FUN[[1]][[1]]@structure[[th.pos1[[1]][[1]][j]]]@dimension, nrow = n)

        fbarre <- stringr::str_replace_all(fbarre, "z", "th")

        res[[FUN[[1]][[1]]@structure[[th.pos1[[1]][[1]][j]]]@arg[1]]] <- eval(parse(text = fbarre))
      }

      if (length(m.pos1[[1]][[1]]) != 0){
        for (j in 1:length(m.pos1[[1]][[1]]))
        {
          FUN[[1]][[1]]@PGF <- stringr::str_replace_all(FUN[[1]][[1]]@PGF, FUN[[1]][[1]]@Param, FUN[[1]][[1]]@parameter)

          param <- paste(FUN[[1]][[1]]@structure[[m.pos1[[1]][[1]][j]]]@Param, num1[[1]][m.pos1[[1]][[1]][j]], sep = "")

          FUN[[1]][[1]]@structure[[m.pos1[[1]][[1]][j]]]@PGF <- stringr::str_replace_all(FUN[[1]][[1]]@structure[[m.pos1[[1]][[1]][j]]]@PGF, FUN[[1]][[1]]@structure[[m.pos1[[1]][[1]][j]]]@Param, param)

          FUN[[1]][[1]]@structure[[m.pos1[[1]][[1]][j]]]@PGF <- stringr::str_replace_all(FUN[[1]][[1]]@PGF, "z", FUN[[1]][[1]]@structure[[m.pos1[[1]][[1]][j]]]@PGF)

          FUN[[1]][[1]]@structure[[m.pos1[[1]][[1]][j]]]@theta <- vapply(FUN[[1]][[1]]@theta, function(x) sum(FUN[[1]][[1]]@structure[[m.pos1[[1]][[1]][j]]]@simul(x, FUN[[1]][[1]]@structure[[m.pos1[[1]][[1]][j]]]@parameter)), rep(0, 1))

          FUN[[1]][[1]]@structure[[m.pos1[[1]][[1]][j]]]@Param <- param

          FUN[[i + 1]][[j]] <- FUN[[1]][[1]]@structure[[m.pos1[[1]][[1]][j]]]
        }}

    }

    else
    {
      th.pos1[[i]] <- list()
      m.pos1[[i]] <- list()

      for (k in 1:length(m.pos1[[i - 1]][[1]]))
      {
        num1[[i]] <- list()
        num1[[i]][[k]] <- paste(i, k, 1:FUN[[i]][[k]]@dimension, sep = "")

        typ <- numeric(FUN[[i]][[k]]@dimension)
        for (j in 1:length(FUN[[i]][[k]]@structure))
          typ[j] <- FUN[[i]][[k]]@structure[[j]]@type

        th.pos1[[i]][[k]] <- which(typ == "Child")
        m.pos1[[i]][[k]] <- which(typ == "Mother")

        if (FUN[[i]][[k]]@arg != 0)
        {
          for (j in 1:length(FUN[[1]][[1]]@arg))
          {
            laplace <- stringr::str_replace_all(FUN[[i]][[k]]@Laplace, FUN[[i]][[k]]@Param, FUN[[i]][[k]]@parameter)
            th <- "-log(runif(n)) / FUN[[i]][[k]]@theta"
            laplace <- parse(text = stringr::str_replace_all(laplace, "z", th))
            res[[FUN[[i]][[k]]@arg[j]]] <- eval(laplace)
          }
        }

        for (j in 1:length(th.pos1[[i]][[k]]))
        {
          param <- paste(FUN[[i]][[k]]@structure[[th.pos1[[i]][[k]][j]]]@Param, num1[[i]][[k]][th.pos1[[i]][[k]][j]], sep = "")
          laplace <- stringr::str_replace_all(FUN[[i]][[k]]@structure[[th.pos1[[i]][[k]][j]]]@Laplace, FUN[[i]][[k]]@structure[[th.pos1[[i]][[k]][j]]]@Param, param)
          fbarre <- stringr::str_replace_all(FUN[[i]][[k]]@PGF, "z", laplace)

          fbarre <- stringr::str_replace_all(fbarre, param, FUN[[i]][[k]]@structure[[th.pos1[[i]][[k]][j]]]@parameter)

          fbarre <- stringr::str_replace_all(fbarre, FUN[[i]][[k]]@Param, FUN[[i]][[k]]@parameter)

          th <- matrix(-log(runif(FUN[[i]][[k]]@structure[[th.pos1[[i]][[k]][j]]]@dimension * n)) / vapply(FUN[[i]][[k]]@theta, function(x) sum(FUN[[i]][[k]]@structure[[th.pos1[[i]][[k]][j]]]@simul(x, FUN[[i]][[k]]@structure[[th.pos1[[i]][[k]][j]]]@parameter)), rep(0, 1)), ncol = FUN[[i]][[k]]@structure[[th.pos1[[i]][[k]][j]]]@dimension, nrow = n)

          fbarre <- stringr::str_replace_all(fbarre, "z", "th")

          res[[FUN[[i]][[k]]@structure[[th.pos1[[i]][[1]][j]]]@arg[1]]] <- eval(parse(text = fbarre))

        }

        if (length(m.pos1[[i]][[k]]) != 0){
          for (j in 1:length(m.pos1[[i]][[k]]))
          {
            FUN[[i]][[k]]@PGF <- stringr::str_replace_all(FUN[[i]][[k]]@PGF, FUN[[i]][[k]]@Param, FUN[[i]][[k]]@parameter)

            param <- paste(FUN[[i]][[k]]@structure[[m.pos1[[i]][[k]][j]]]@Param, num1[[i]][[k]][m.pos1[[i]][[k]][j]], sep = "")

            FUN[[i]][[k]]@structure[[m.pos1[[i]][[k]][j]]]@PGF <- stringr::str_replace_all(FUN[[i]][[k]]@structure[[m.pos1[[i]][[k]][j]]]@PGF, FUN[[i]][[k]]@structure[[m.pos1[[i]][[k]][j]]]@Param, param)

            FUN[[i]][[k]]@structure[[m.pos1[[i]][[k]][j]]]@PGF <- stringr::str_replace_all(FUN[[i]][[k]]@PGF, "z", FUN[[i]][[k]]@structure[[m.pos1[[i]][[k]][j]]]@PGF)

            FUN[[i]][[k]]@structure[[m.pos1[[i]][[k]][j]]]@theta <- vapply(FUN[[i]][[k]]@theta, function(x) sum(FUN[[i]][[k]]@structure[[m.pos1[[i]][[k]][j]]]@simul(x, FUN[[i]][[k]]@structure[[m.pos1[[i]][[k]][j]]]@parameter)), rep(0, 1))

            FUN[[i]][[k]]@structure[[m.pos1[[i]][[k]][j]]]@Param <- param

            FUN[[i + 1]][[j]] <- FUN[[i]][[k]]@structure[[m.pos1[[i]][[k]][j]]]

          }}
      }
    }

  }
  do.call(cbind, res)
})

#' @rdname dCompCop
#' @export

rCompCop2 <- compiler::cmpfun(function(n, FUN, level){

  if (FUN@type == "Child")
  {
    FUN@theta <- FUN@simul(n, FUN@parameter)
    laplace <- stringr::str_replace_all(FUN@Laplace, FUN@Param, FUN@parameter)
    th <- "-log(runif(FUN@dimension * n)) / FUN@theta"
    laplace <- parse(text = stringr::str_replace_all(laplace, "z", th))
    return(matrix(eval(laplace), ncol = FUN@dimension, nrow = n))
  }

  t <- FUN
  res <- list()
  FUN <- list()

  th.pos1 <- list()
  m.pos1 <- list()

  th.pos1[[1]] <- list()
  m.pos1[[1]] <- list()

  FUN[[1]] <- list()
  FUN[[1]][[1]] <- t


  for (i in 1:level)
  {
    FUN[[i + 1]] <- list()

    if (i == 1){

      FUN[[1]][[1]]@theta <- FUN[[1]][[1]]@simul(n, FUN[[1]][[1]]@parameter)

      FUN[[1]][[1]]@PGF <- stringr::str_replace_all(FUN[[1]][[1]]@PGF, FUN[[1]][[1]]@Param, FUN[[1]][[1]]@parameter)

      if (FUN[[1]][[1]]@type == "Mother"){
      typ <- numeric(FUN[[1]][[1]]@dimension)
      for (j in 1:length(FUN[[1]][[1]]@structure))
        typ[j] <- FUN[[1]][[1]]@structure[[j]]@type}

      th.pos1[[1]][[1]] <- which(typ == "Child")
      m.pos1[[1]][[1]] <- which(typ == "Mother")

      if (length(FUN[[1]][[1]]@arg) > 1 || {length(FUN[[1]][[1]]@arg) == 1 && FUN[[1]][[1]]@arg != 0})
      {
        for (j in 1:length(FUN[[1]][[1]]@arg))
        {
          laplace <- stringr::str_replace_all(FUN[[1]][[1]]@Laplace, FUN[[1]][[1]]@Param, FUN[[1]][[1]]@parameter)
          th <- "-log(runif(n)) / FUN[[1]][[1]]@theta"
          laplace <- parse(text = stringr::str_replace_all(laplace, "z", th))
          res[[FUN[[1]][[1]]@arg[j]]] <- eval(laplace)
        }
      }

      if (length(th.pos1[[1]][[1]]) != 0){
      for (j in 1:length(th.pos1[[1]][[1]]))
      {
        laplace <- stringr::str_replace_all(FUN[[1]][[1]]@structure[[th.pos1[[1]][[1]][j]]]@Laplace, FUN[[1]][[1]]@structure[[th.pos1[[1]][[1]][j]]]@Param, FUN[[1]][[1]]@structure[[th.pos1[[1]][[1]][j]]]@parameter)
        fbarre <- stringr::str_replace_all(FUN[[1]][[1]]@PGF, "z", laplace)

        argg <- FUN[[1]][[1]]@structure[[th.pos1[[1]][[1]][j]]]@arg

        th2 <- vapply(FUN[[1]][[1]]@theta, function(x) sum(FUN[[1]][[1]]@structure[[th.pos1[[1]][[1]][j]]]@simul(x, FUN[[1]][[1]]@structure[[th.pos1[[1]][[1]][j]]]@parameter)), rep(0, 1))

        res2 <- list()
        for (l in 1:length(argg))
        {
          th <- -log(runif(n)) / th2

          fbarre <- stringr::str_replace_all(fbarre, "z", "th")

          res[[argg[l]]] <- eval(parse(text = fbarre))
        }
      }}

      if (length(m.pos1[[1]][[1]]) != 0){
        for (j in 1:length(m.pos1[[1]][[1]]))
        {

          FUN[[1]][[1]]@structure[[m.pos1[[1]][[1]][j]]]@PGF <- stringr::str_replace_all(FUN[[1]][[1]]@structure[[m.pos1[[1]][[1]][j]]]@PGF, FUN[[1]][[1]]@structure[[m.pos1[[1]][[1]][j]]]@Param, FUN[[1]][[1]]@structure[[m.pos1[[1]][[1]][j]]]@parameter)

          FUN[[1]][[1]]@structure[[m.pos1[[1]][[1]][j]]]@PGF <- stringr::str_replace_all(FUN[[1]][[1]]@PGF, "z", FUN[[1]][[1]]@structure[[m.pos1[[1]][[1]][j]]]@PGF)

          FUN[[1]][[1]]@structure[[m.pos1[[1]][[1]][j]]]@theta <- vapply(FUN[[1]][[1]]@theta, function(x) sum(FUN[[1]][[1]]@structure[[m.pos1[[1]][[1]][j]]]@simul(x, FUN[[1]][[1]]@structure[[m.pos1[[1]][[1]][j]]]@parameter)), rep(0, 1))

          FUN[[i + 1]][[j]] <- FUN[[1]][[1]]@structure[[m.pos1[[1]][[1]][j]]]
        }}

    }

    else
    {
      th.pos1[[i]] <- list()
      m.pos1[[i]] <- list()

      for (k in 1:length(m.pos1[[i - 1]][[1]]))
      {

        typ <- numeric(FUN[[i]][[k]]@dimension)
        for (j in 1:length(FUN[[i]][[k]]@structure))
          typ[j] <- FUN[[i]][[k]]@structure[[j]]@type

        th.pos1[[i]][[k]] <- which(typ == "Child")
        m.pos1[[i]][[k]] <- which(typ == "Mother")

        if (FUN[[i]][[k]]@arg != 0)
        {
          for (j in 1:length(FUN[[i]][[k]]@arg))
          {
            laplace <- stringr::str_replace_all(FUN[[i]][[k]]@Laplace, FUN[[i]][[k]]@Param, FUN[[i]][[k]]@parameter)
            th <- "-log(runif(n)) / FUN[[i]][[k]]@theta"
            laplace <- parse(text = stringr::str_replace_all(laplace, "z", th))
            res[[FUN[[i]][[k]]@arg[j]]] <- eval(laplace)
          }
        }

        if (length(th.pos1[[i]][[k]]) != 0){
        for (j in 1:length(th.pos1[[i]][[k]]))
        {
          laplace <- stringr::str_replace_all(FUN[[i]][[k]]@structure[[th.pos1[[i]][[k]][j]]]@Laplace, FUN[[i]][[k]]@structure[[th.pos1[[i]][[k]][j]]]@Param, FUN[[i]][[k]]@structure[[th.pos1[[i]][[k]][j]]]@parameter)

          fbarre <- stringr::str_replace_all(FUN[[i]][[k]]@PGF, "z", laplace)

          argg <- FUN[[i]][[k]]@structure[[th.pos1[[i]][[1]][j]]]@arg

          th2 <- vapply(FUN[[i]][[k]]@theta, function(x) sum(FUN[[i]][[k]]@structure[[th.pos1[[i]][[k]][j]]]@simul(x, FUN[[i]][[k]]@structure[[th.pos1[[i]][[k]][j]]]@parameter)), rep(0, 1))

          res2 <- list()
          for (l in 1:length(argg))
          {
            th <- -log(runif(n)) / th2

            fbarre <- parse(text = stringr::str_replace_all(fbarre, "z", "th"))

            res[[argg[l]]] <- eval(fbarre)
          }

          #res[[argg[1]]] <- matrix(unlist(res2), ncol = FUN[[1]][[1]]@structure[[th.pos1[[1]][[1]][j]]]@dimension, nrow = n)
        }}

        if (length(m.pos1[[i]][[k]]) != 0){
          for (j in 1:length(m.pos1[[i]][[k]]))
          {

            FUN[[i]][[k]]@structure[[m.pos1[[i]][[k]][j]]]@PGF <- stringr::str_replace_all(FUN[[i]][[k]]@structure[[m.pos1[[i]][[k]][j]]]@PGF,
                                                                                           FUN[[i]][[k]]@structure[[m.pos1[[i]][[k]][j]]]@Param,
                                                                                           FUN[[i]][[k]]@structure[[m.pos1[[i]][[k]][j]]]@parameter)

            FUN[[i]][[k]]@structure[[m.pos1[[i]][[k]][j]]]@PGF <- stringr::str_replace_all(FUN[[i]][[k]]@PGF, "z", FUN[[i]][[k]]@structure[[m.pos1[[i]][[k]][j]]]@PGF)

            FUN[[i]][[k]]@structure[[m.pos1[[i]][[k]][j]]]@theta <- vapply(FUN[[i]][[k]]@theta, function(x) sum(FUN[[i]][[k]]@structure[[m.pos1[[i]][[k]][j]]]@simul(x, FUN[[i]][[k]]@structure[[m.pos1[[i]][[k]][j]]]@parameter)), rep(0, 1))

            FUN[[i + 1]][[j]] <- FUN[[i]][[k]]@structure[[m.pos1[[i]][[k]][j]]]

          }}
      }
    }

  }
  do.call(cbind, res)
})

#' Bounded Probability for Random Vectors
#'
#' Probability computation
#' @param FUN Multivariate CDF/Survival function
#' @param dim Dimension of the random vector
#' @param data Matrix of bounds (by row)
#' @param prob FUN is a 'CDF' or 'Survival' function ? (either 'CDF' or 'Survival')
#' @return The probability
#' @importFrom mgcv uniquecombs
#' @importFrom gtools permutations
#' @export

f.mult <- compiler::cmpfun(function(FUN, dim, data, prob = "CDF")
{
  d <- dim
  L <- list()

  for(i in 0:d)
  {
    L[[i+1]] <- mgcv::uniquecombs(gtools::permutations(d, d, set = FALSE, v = c(rep(1, i), rep(2, d - i))))
  }
  pos <- do.call(rbind, L)

  sign <- list()
  for (y in 0:d)
    sign[[y + 1]] <- (-1)^((y + 1) %% 2 == 0) * rep(1, choose(d, y))

  if (prob == "Survival")
    sign <- unlist(sign) * (-1)^(dim %% 2 != 0) ## À VÉRIFIER
  if (prob == "CDF")
    sign <- unlist(sign)
  else
    stop("Wrong 'prob' input")

  sub1 <- numeric(length(sign))
  for (z in 1:length(sign))
  {
    position <- pos[z,]

    sub2 <- numeric(d)
    for (k in 1:d)
      sub2[k] <- data[k,][position[k]]
    sub1[z] <- FUN(sub2)
  }

  sum(sub1 * sign)
})

#' Estimation for one level hierarchical copulas
#'
#'
#' @param data Data used for the estimation
#' @param struc Known structure (S4)
#' @param lower Lower bound vector for optimization (from left to right, per level)
#' @param upper Upper bound vector for optimization (from left to right, per level)
#' @param der List of derivatives (from left to right, per level)
#' @export

CompCopEstim <- compiler::cmpfun(function(data, struc, lower, upper, der)
{
  res <- data
  kk <- 1 ## Compteur pour les bornes

  ## lower: Borne inf. pour l'optimization -> de gauche à droite, par niveau
  ## upper: Borne sup. pour l'optimization -> de gauche à droite, par niveau

  ## 1: Sortir le nombre d'éléments dans chaque sous-groupe
  n <- numeric(length(struc@structure))
  for (i in 1:length(struc@structure))
    n[i] <- length(struc@structure[[i]]@arg)

  uu <- paste("x", 1:length(n), sep = "")

  MAT1 <- list()
  for (i in 1:length(n))
  {
    for (b in 1:length(struc@structure[[i]]@arg))
      MAT1[[i]] <- struc@structure[[i]]@arg
  }

  ini <- paste("1:n[", 1:length(n), "]", sep = "", collapse = ", ")
  out <- "expand.grid(z)"
  MAT <- eval(parse(text = stringr::str_replace_all(out, "z", ini)))

  grp <- length(n)
  xx1 <- paste("x", 1:length(n), sep = "")
  xx <- paste(xx1, collapse = ", ")
  test <- "function(z)"
  test <- paste(stringr::str_replace_all(test, "z", xx), " c(z)", sep = "")
  xx2 <- paste("MAT1[[", 1:length(n), "]][", xx1, "]", sep = "")
  xx2 <- paste(xx2, collapse = ", ")
  test <- stringr::str_replace_all(test, "z", xx2)
  test <- parse(text = test)

  input <- paste("MAT[,", 1:length(MAT[1,]), "]", collapse = ", ", sep = "")
  final <- "mapply(eval(test), z)"
  final <- stringr::str_replace_all(final, "z", input)
  MAT <- t(eval(parse(text = final)))

  ## 4: Optimiser pour le param de M
  l <- list()
  for (i in 1:length(MAT[,1]))
    l[[i]] <- res[,MAT[i,]]

  logvM <- function(par, data)
  {
    alpha <- par
    res <- rep(1, length(data[[1]][,1]))
    #ini <- paste("data[[i]][,", 1:length(n), "]", sep = "", collapse = ", ")
    for (i in 1:length(data))
    {
      for (b in 1:length(uu))
        eval(parse(text = paste("x", b, " <- data[[i]][,", b, "]", sep = "")))
      res <- res * eval(der[[1]])
    }
    -sum(log(res))
  }

  alpha0 <- optimize(logvM, c(lower[kk], upper[kk]), data = l)$minimum
  kk <- kk + 1

  ## 5: Estimer les paramètres des sous-groupes
  test <- struc
  gamma <- numeric(length(n))
  for (j in 1:length(n))
  {
    logvB <- function(par, data)
    {
      gamma <- alpha0
      alpha <- par
      for (b in 1:length(test@structure[[j]]@arg))
        eval(parse(text = paste("u", b, " <- data[,", b, "]", sep = "")))
      -sum(log(eval(der[[j + 1]])))
    }

    gamma[j] <- suppressWarnings(optimize(logvB, c(lower[kk], upper[kk]), data = res[,MAT1[[j]]])$minimum)
    kk <- kk + 1
  }
  list("M" = alpha0,
       "B" = gamma)
})

#' Estimation and bootstrap for one level hierarchical copulas
#'
#'
#' @param m Number of bootstraps
#' @param nn size of the sample per bootstrap
#' @param struc Known structure (S4)
#' @param lower Lower bound vector for optimization (from left to right, per level)
#' @param upper Upper bound vector for optimization (from left to right, per level)
#' @param data Data used for the estimation if their is no bootstrap
#' @export

CompCopBootstrap <- function(m, nn, struc, lower, upper, data = NULL)
{
  vk <- 1:length(struc@structure)
  uB <- paste("B", vk, sep = "")
  ll <- list()

  ## M ##
  expr <- pCop(struc@cop(0.99, length(vk)), vector = F, express = T) ## 0.99 a aucun impact
  uu <- paste("x", 1:length(uB), sep = "")
  for (i in 1:length(uB))
    expr <- Deriv::Deriv(expr, uu[i], cache.exp = F)
  ll[[1]] <- expr

  ## B ##
  for (j in 1:length(uB))
  {
    Phi <- stringr::str_replace_all(struc@PGF, "z", struc@structure[[j]]@Laplace)
    yy <- paste("y", 1:length(struc@structure[[j]]@arg), sep = "")
    uu <- paste("u", 1:length(struc@structure[[j]]@arg), sep = "")
    ini <- paste(yy, collapse = " + ")
    PhiInv <- stringr::str_replace_all(struc@structure[[j]]@LaplaceInv, "z", struc@PGFInv)
    for (b in 1:length(struc@structure[[j]]@arg))
    {
      ini <- stringr::str_replace_all(ini, yy[b], PhiInv)
      ini <- stringr::str_replace_all(ini, "z", uu[b])
    }

    cop <- parse(text = stringr::str_replace_all(Phi, "z", ini))
    tt <- eval(parse(text = paste("{for (b in 1:length(struc@structure[[", j, "]]@arg))
                                  cop <- Deriv::Deriv(cop, uu[b], cache.exp = F)}", sep = "")))
    ll[[j + 1]] <- cop
  }

  if (is.null(data) == FALSE)
    CompCopEstim(data, struc, lower, upper, ll)
  else
  {
    M <- numeric(m)
    B <- matrix(NA, ncol = length(struc@structure), nrow = m)
    for (i in 1:m)
    {
      x <- CompCopEstim(rCompCop2(nn, struc, 1), struc, lower, upper, der = ll)
      M[i] <- x$M
      B[i,] <- x[[2]]
    }
    list("M" = M,
         "B" = B)
  }
}

#' Obtain a node with its genetic code
#'
#' @param path Genetic code of the node
#' @param str Structure
#'
#' @export

Node <- function(path, str)
{
  if (length(path) == 1)
    str
  else
  {
    if (path[2] == 0)
      stop("There is no node associated to 0")

    struc <- str@structure[[path[2]]]
    Node(path[-2], struc)
  }
}

#' Derivative of the power of a pgf
#'
#' @param tt What should be the expression evaluated to ? (the variable)
#' @param n The number of derivatives
#' @param s The power
#' @param str The structure
#'
#' @export

PGF.Power <- function(tt, n, s, str)
{
  ## Trouver les combinaisons
  comb <- "test <- expand.grid(z)"
  ini <- numeric(s)
  for (i in 1:s)
    ini[i] <- "0:n"
  ini <- paste(ini, collapse = ", ")
  eval(parse(text = stringr::str_replace_all(comb, "z", ini)))

  test <- as.matrix(test)
  test2 <- matrix(NA, ncol = s, nrow = length(test[,1]))
  kk <- 1
  for (i in 1:length(test[,1]))
  {
    if (sum(test[i,]) == n)
    {
      test2[kk,] <- test[i,]
      kk <- kk + 1
    }
  }
  test2 <- test2[-(kk:length(test[,1])),]
  test2 <- as.matrix(test2)

  ## Construire la dérivée
  ini <- numeric(length(test2[,1]))
  xx2 <- paste("x", s:1, sep = "")
  for (i in 1:length(test2[,1]))
  {
    xx1 <- paste("(x", s:1, ")", sep = "", collapse = " * ")
    for (j in 1:s)
    {
      xx1 <- stringr::str_replace_all(xx1, xx2[j], str@Der(tt, test2[i, j], "PGF"))
    }
    res1 <- factorial(n) / prod(sapply(1:s, function(t) factorial(test2[i, t])))
    res2 <- xx1

    ini[i] <- paste("(", res1, ") * (", res2, ")", sep = "")
  }

  paste("(", ini, ")", sep = "", collapse = " + ")
}

#' Expression for the density of a Theta (level 1)
#'
#' @param str Structure
#' @param grp Position of the group
#' @param dim Dimension of the group
#' @param simplify Should the expression be simply ? (may take time)
#'
#' @details The goal is to soon introduce genetic codes to simplify the process
#'
#' @export

density_level1_EXPR <- function(str, grp, dim, simplify = FALSE)
{
  uu <- paste("u", 1:dim, sep = "")
  yy <- paste("y", 1:dim, sep = "")
  expr1 <- paste("y", 1:dim, sep = "", collapse = " * ")
  nu <- paste(yy, collapse = " + ")
  for (i in 1:dim)
  {
    ini <- paste("(", stringr::str_replace_all(str@structure[[grp]]@Der(str@PGFInv, 1, "LaplaceInv"), "z", uu[i]), ")", sep = "")
    ini <- paste(ini, " * (", str@Der(uu[i], 1, "PGFInv"), ")", sep = "")
    expr1 <- stringr::str_replace_all(expr1, yy[i], ini)

    ini <- stringr::str_replace_all(str@structure[[grp]]@LaplaceInv, "z", str@PGFInv)
    nu <- stringr::str_replace_all(nu, yy[i], ini)
    nu <- stringr::str_replace_all(nu, "z", uu[i])
  }

  res1 <- numeric(dim)
  for (r in 1:dim)
  {
    expr2 <- str@Der(stringr::str_replace_all(str@structure[[grp]]@Laplace, "z", nu), r, "PGF")
    res2 <- numeric(r)
    for (s in 1:r)
    {
      res2[s] <- paste("((-1)^(", r - s, ") / (", factorial(s) * factorial(r - s), ") * (",
                       stringr::str_replace_all(str@structure[[grp]]@Laplace, "z", nu), ")^(", r - s, ") * (",
                       stringr::str_replace_all(stringr::str_replace_all(str@structure[[grp]]@Der("z", dim, "Laplace"), "alpha", paste("(alpha * ", s, ")", sep = "")),
                                                "z", nu), "))", sep = "")
    }
    res1[r] <- paste("(", expr2, ") * (", paste(res2, collapse = " + "), ")", sep = "")
  }
  res <- paste(res1, collapse = " + ")

  expr.final <- paste("(", expr1, ") * (", res, ")", sep = "")

  if (simplify == TRUE)
    Deriv::Simplify(expr.final)
  else
    expr.final
}

#' Expression for the density of a Theta (level 2)
#'
#' @param str Structure
#' @param grpM Position of the M
#' @param grpB Position of the Theta (under M)
#' @param dim Dimension of the group
#' @param simplify Should the expression be simply ? (may take time)
#'
#' @details The goal is to soon introduce genetic codes to simplify the process
#'
#' @export

density_level2_EXPR <- function(str, grpM, grpB, dim, simplify = FALSE)
{
  dimm <- dim
  M0 <- str
  M1 <- M0@structure[[grpM]]
  B11 <- M1@structure[[grpB]]

  ini1 <- B11@Der(stringr::str_replace_all(M1@PGFInv, "gamma", "gamma1"),
                  1, "LaplaceInv")
  ini1 <- stringr::str_replace_all(ini1, "z",
                                   stringr::str_replace_all(M0@PGFInv, "gamma", "gamma0"))
  ini2 <- stringr::str_replace_all(M1@Der("z", 1, "PGFInv"), "gamma", "gamma1")
  ini2 <- stringr::str_replace_all(ini2, "z",
                                   stringr::str_replace_all(M0@PGFInv, "gamma", "gamma0"))
  ini3 <- stringr::str_replace_all(M0@Der("z", 1, "PGFInv"), "gamma", "gamma0")
  ini <- paste("(", ini1, ") * (", ini2, ") * (", ini3, ")", sep = "")

  uu <- paste("u", 1:dimm, sep = "")

  expr1 <- numeric(dimm)
  for (i in 1:dimm)
    expr1[i] <- stringr::str_replace_all(ini, "z", uu[i])

  expr1 <- paste("(", expr1, ")", sep = "", collapse = " * ")
  expr1 <- Deriv::Simplify(expr1)

  ## BLEU ET ROUGE ##

  nu <- stringr::str_replace_all(B11@LaplaceInv, "z",
                                 stringr::str_replace_all(M1@PGFInv, "gamma", "gamma1"))
  nu <- stringr::str_replace_all(nu, "z",
                                 stringr::str_replace_all(M0@PGFInv, "gamma", "gamma0"))
  nuu <- numeric(dimm)
  for (i in 1:dimm)
    nuu[i] <- stringr::str_replace_all(nu, "z", uu[i])
  nu <- paste("(", nuu, ")", sep = "", collapse = " + ")
  nu <- Deriv::Simplify(nu)

  input1 <- stringr::str_replace_all(M1@PGF, "gamma", "gamma1")
  input1 <- stringr::str_replace_all(input1, "z", B11@Laplace)
  input1 <- stringr::str_replace_all(input1, "z", nu)

  vec1 <- numeric(dimm)
  for (r in 1:dimm)
  {
    vec21 <- numeric(r)
    for (i in 1:r)
    {
      ini <- stringr::str_replace_all(M0@Der("z", i, "PGF"),
                                      "gamma", "gamma0")
      ini <- stringr::str_replace_all(ini, "z", input1)

      vec3 <- numeric(i)
      for (j in 1:i)
      {
        coeff1 <- (-1)^(i - j) / factorial(j) / factorial(i - j)

        coeff2 <- input1
        coeff2 <- paste("(", coeff2, ")^(", i - j, ")", sep = "")

        coeff3 <- PGF.Power("z", r, j, M1)
        coeff3 <- stringr::str_replace_all(coeff3, "gamma", "gamma1")
        coeff3 <- stringr::str_replace_all(coeff3, "z",
                                           stringr::str_replace_all(B11@Laplace, "z", nu))
        coeff <- paste("(", coeff1, ") * (", coeff2, ") * (", coeff3, ")", sep = "")

        vec3[j] <- coeff
      }

      vec21[i] <- paste("(", ini, ") * (",
                        paste("(", vec3, ")", sep = "", collapse = " + "), ")", sep = "")
    }

    allo1 <- paste("(", vec21, ")", sep = "", collapse = " + ")

    vec22 <- numeric(r)
    for (s in 1:r)
    {
      coeff1 <- (-1)^(r - s) / factorial(s) / factorial(r - s)

      coeff2 <- paste("(", stringr::str_replace_all(B11@Laplace, "z", nu), ")^(",
                      r - s, ")", sep = "")

      coeff3 <- B11@Der("z", dimm, "Laplace")
      coeff3 <- stringr::str_replace_all(coeff3, "alpha", paste("(alpha", " * ", s, ")", sep = ""))
      coeff3 <- stringr::str_replace_all(coeff3, "z", nu)

      coeff <- paste("(", coeff1, ") * (", coeff2, ") * (", coeff3, ")", sep = "")

      vec22[s] <- coeff
    }

    allo2 <- paste("(", vec22, ")", sep = "", collapse = " + ")
    vec1[r] <- paste("(", allo1, ") * (", allo2, ")", sep = "")
  }

  final <- paste("(", vec1, ")", sep = "", collapse = " + ")
  final <- paste("(", expr1, ") * (", final, ")", sep = "")

  if (simplify)
    Deriv::Simplify(final)
  else
    final
}
