#' LST of a Child Node
#'
#' With a specific path and a predefined structure (S4 class of a type 'Mother'), t
#' his function returns the LST expression of the corresponding node with a
#' specific variable.
#'
#' @description The function Lap_Child() ...
#' @param code Genetic code of the child node (can be a leaf i.e. end by 0)
#' @param str Object of class Mother (the structure)
#' @param tt Output variable to be used ('z' by default)
#'
#' @rdname Lap_Child
#'
#' @seealso \link{InvLap_Child}
#'
#'
#' @examples
#'
#' str <- GEO(0.1, NULL, list(GAMMA(0.1, 1:2, NULL),
#'                            GAMMA(0.2, 3:4, NULL)))
#'
#' InvLap_Child(c(0,2), str)
#'
#' @author Simon-Pierre Gadoury
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
