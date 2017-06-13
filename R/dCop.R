#' Density for Archimedean Copulas Objects
#'
#' @param copula An Archimedean copula class object
#' @param vector If false, returns a function with (x_1, x_2, ..., x_dim, alpha) as arguments.
#' @param express If true, returns an expression.
#' @param code If true, copies the LaTeX code to clipboard.
#' @param operator Type of cumputer used (only necessary in the case of internal problem)
#' @return Either an expression, function or code.
#' @author Simon-Pierre Gadoury
#'
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
