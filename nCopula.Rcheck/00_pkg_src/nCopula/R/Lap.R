#' LST of a Node
#'
#' @description With a specific path and a predefined structure
#' (S4 class of a type 'Mother'), returns the Laplace-Stieltjes Transform expression of
#' the corresponding node with a specific variable.
#'
#' @param code genetic code (numeric vector) of the node (can be a leaf i.e. end by 0).
#' @param structure object of class Mother (the structure).
#' @param outVar output variable to be used ('z' by default).
#' @param par Should the parameters be values ('value') or variables ('variable') ?
#'
#' @rdname Lap
#'
#' @seealso \link{InvLap}
#'
#' @importFrom utils head
#'
#' @return A character string giving the LST of the specified node.
#'
#' @details For mother nodes, parameters are always called 'gamma' and for child nodes, parameters are
#' always called 'alpha'. Furthermore, to recognize the parameters, the path is inserted at the end.
#' For exemple, a child node with path (0,2,1) will have the parameter 'alpha021'.
#'
#' @examples
#'
#' structure <- GEO(0.1, NULL, list(GAMMA(0.1, 1:2, NULL),
#'                            GAMMA(0.2, 3:4, NULL)))
#'
#' Lap(c(0,2), structure, outVar = 'z', par = 'value')
#'
#' @author Simon-Pierre Gadoury
#' @export

Lap <- function(code, structure, outVar = "z", par = "value")
{
  str_ini <- structure

  if (par == "value")
  {
    lap <- stringr::str_replace_all(structure@PGF, structure@Param, structure@parameter)
    for (i in 2:length(code))
    {
      code2 <- head(code, i)
      str2 <- Node(code2, str_ini)

      if (str2@type == "Mother" && i != length(code))
        ini <- stringr::str_replace_all(str2@PGF, str2@Param, str2@parameter)
      else
        ini <- stringr::str_replace_all(str2@Laplace, str2@Param, str2@parameter)

      lap <- stringr::str_replace_all(lap, "z", ini)
    }
  }
  else if (par == "variable")
  {
    lap <- stringr::str_replace_all(structure@PGF, structure@Param, paste(structure@Param, "0", sep = ""))
    for (i in 2:length(code))
    {
      code2 <- head(code, i)
      str2 <- Node(code2, str_ini)

      if (str2@type == "Mother" && i != length(code))
        ini <- stringr::str_replace_all(str2@PGF, str2@Param,
                                        paste(str2@Param, paste(code2, collapse = ""), sep = ""))
      else
        ini <- stringr::str_replace_all(str2@Laplace, str2@Param,
                                        paste(str2@Param, paste(code2, collapse = ""), sep = ""))

      lap <- stringr::str_replace_all(lap, "z", ini)
    }
  }
  stringr::str_replace_all(lap, "z", outVar)
}
