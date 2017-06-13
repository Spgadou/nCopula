#' Density, Cdf, and Random Number Generator for Copulas Constructed Through Compounding
#'
#' @param FUN Object of class Mother
#' @param func If true, returns a function
#' @param code If true, copies the LaTeX code to the clipboard
#' @param operator Type of cumputer used (only necessary in the case of internal problem)
#'
#' @details rCompCop2 is more general (and easier to use) than rCompCop, but is slower.
#'
#' @author Simon-Pierre Gadoury
#'
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
