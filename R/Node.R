#' Obtain a node with its genetic code
#'
#' @param path Genetic code of the node
#' @param str The structure
#'
#' @author Simon-Pierre Gadoury
#'
#' @export

Node <- function(path, str)
{
  if (length(path) == 1)
    str
  else
  {
    if (path[2] == 0)
      str
    else
    {
      struc <- str@structure[[path[2]]]
      Node(path[-2], struc)
    }
  }
}
