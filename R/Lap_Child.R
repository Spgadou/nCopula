#' LST of a Child Node
#'
#' @param code Genetic code of the child node (can be a leaf i.e. end by 0)
#' @param str Object of class Mother (the structure)
#' @param tt Output variable to be used ('z' by default)
#'
#' @author Simon-Pierre Gadoury
#'
#' @export

Lap_Child <- function(code, str, tt = "z")
{
  str_ini <- str
  lap <- stringr::str_replace_all(str@PGF, str@Param, str@parameter)
  for (i in 2:length(code))
  {
    code2 <- head(code, i)
    str2 <- Node(code2, str_ini)

    if (str2@type == "Mother")
      ini <- stringr::str_replace_all(str2@PGF, str2@Param, str2@parameter)
    else
      ini <- stringr::str_replace_all(str2@Laplace, str2@Param, str2@parameter)

    lap <- stringr::str_replace_all(lap, "z", ini)
  }
  stringr::str_replace_all(lap, "z", tt)
}
