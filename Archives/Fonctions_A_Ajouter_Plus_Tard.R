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
