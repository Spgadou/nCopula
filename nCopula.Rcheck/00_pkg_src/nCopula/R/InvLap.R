#' Inverse LST of a Node
#'
#' @description With a specific path and a predefined structure
#' (S4 class of a type 'Mother'), returns the inverse Laplace-Stieltjes Transform expression of
#' the corresponding node with a specific variable.
#'
#' @param code the genetic code (numeric vector) of the node (can be a leaf i.e. end by 0).
#' @param structure an object of class Mother (the structure).
#' @param outVar the output variable to be used ('z' by default).
#' @param par logical. Should the parameters be values ('value') or variables ('variable') ?
#'
#' @seealso \link{Lap}
#'
#' @importFrom utils tail
#'
#' @rdname InvLap
#'
#' @return A character string giving the inverse LST of the specified node.
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
#' InvLap(c(0,2), structure, outVar = 'z', par = 'value')
#'
#' @author Simon-Pierre Gadoury
#' @export

InvLap <- function(code, structure, outVar = "z", par = "value")
{
  str_next <- Node(code, structure)

  if (par == 'value')
  {
    if (tail(code, 1) != 0)
    {
      ini <- stringr::str_replace_all(str_next@LaplaceInv, str_next@Param, str_next@parameter)
      final <- ini

      for (i in (length(code) - 1):1)
      {
        code2 <- head(code, i)
        str_next <- Node(code2, structure)
        ini <- stringr::str_replace_all(str_next@PGFInv, str_next@Param, str_next@parameter)

        final <- stringr::str_replace_all(final, "z", ini)
      }
    }
    else
    {
      ini <- stringr::str_replace_all(str_next@PGFInv, str_next@Param, str_next@parameter)
      final <- ini

      if (length(code) > 2)
      {
        for (i in (length(code) - 2):1)
        {
          code2 <- head(code, i)
          str_next <- Node(code2, structure)
          ini <- stringr::str_replace_all(str_next@PGFInv, str_next@Param, str_next@parameter)

          final <- stringr::str_replace_all(final, "z", ini)
        }
      }
    }
  }
  else if (par == 'variable')
  {
    if (tail(code, 1) != 0)
    {
      ini <- stringr::str_replace_all(str_next@LaplaceInv, str_next@Param,
                                      paste(str_next@Param, paste(code, collapse = ""), sep = ""))
      final <- ini

      for (i in (length(code) - 1):1)
      {
        code2 <- head(code, i)
        str_next <- Node(code2, structure)
        ini <- stringr::str_replace_all(str_next@PGFInv, str_next@Param,
                                        paste(str_next@Param, paste(code2, collapse = ""), sep = ""))

        final <- stringr::str_replace_all(final, "z", ini)
      }
    }
    else
    {
      ini <- stringr::str_replace_all(str_next@PGFInv, str_next@Param,
                                      paste(str_next@Param, paste(code, collapse = ""), sep = ""))
      final <- ini

      if (length(code) > 2)
      {
        for (i in (length(code) - 2):1)
        {
          code2 <- head(code, i)
          str_next <- Node(code2, structure)
          ini <- stringr::str_replace_all(str_next@PGFInv, str_next@Param,
                                          paste(str_next@Param, paste(code2, collapse = ""), sep = ""))

          final <- stringr::str_replace_all(final, "z", ini)
        }
      }
    }
  }

  final <- stringr::str_replace_all(final, "z", outVar)
  final
}


