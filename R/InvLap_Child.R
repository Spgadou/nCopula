#' Inverse LST of a Child Node
#'
#' @param code Genetic code of the child node
#' @param str Object of class Mother (the structure)
#' @param tt Output variable to be used ('z' by default)
#'
#' @author Simon-Pierre Gadoury
#'
#' @export

InvLap_Child <- function(code, str, tt = "z")
{
  str_next <- Node(code, str)

  if (tail(code, 1) != 0)
  {
    ini <- stringr::str_replace_all(str_next@LaplaceInv, str_next@Param, str_next@parameter)
    final <- ini

    for (i in (length(code) - 1):1)
    {
      code2 <- head(code, i)
      str_next <- Node(code2, str)
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
        str_next <- Node(code2, str)
        ini <- stringr::str_replace_all(str_next@PGFInv, str_next@Param, str_next@parameter)

        final <- stringr::str_replace_all(final, "z", ini)
      }
    }

  }

  final <- stringr::str_replace_all(final, "z", tt)
  final
}


