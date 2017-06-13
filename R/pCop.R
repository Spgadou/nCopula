#' Cdf, and Random Number Generator for Copulas
#'
#' @param copula An Archimedean copula class object
#' @param vector If false, returns a function with (x_1, x_2, ..., x_dim, alpha) as arguments.
#' @param express If true, returns an expression.
#' @param code If true, copies the LaTeX code to clipboard.
#' @param operator Type of cumputer used (only necessary in the case of internal problem)
#' @return Either an expression, function, code, or sampled data.
#'
#' @author Simon-Pierre Gadoury
#'
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
