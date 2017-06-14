#' Obtain a node with its genetic code
#'
#' @description Obtain a node with its genetic code.
#' @param path Genetic code of the node
#' @param str The structure
#'
#' @examples
#' # We directly give the path of the desired node.
#' Node(c(0,2,2), LOG(0.5, NULL, list(GAMMA(1/30, c(5,6), NULL),
#'                               LOG(0.1, NULL, list(GAMMA(1/30, c(1,2), NULL),
#'                               GAMMA(1/30, c(3,4), NULL))))))
#'
#' # Here we provide the path with the GeneticCodes function of this package.
#' Node(GeneticCodes(LOG(0.5, NULL, list(GAMMA(1/30, c(5,6), NULL),
#'                                      LOG(0.1, NULL, list(GAMMA(1/30, c(1,2), NULL),
#'                                      GAMMA(1/30, c(3,4), NULL))))))[[3]],
#'                                  LOG(0.5, NULL, list(GAMMA(1/30, c(5,6), NULL),
#'                                  LOG(0.1, NULL, list(GAMMA(1/30, c(1,2), NULL),
#'                                  GAMMA(1/30, c(3,4), NULL))))))
#'
#' @return The node distribution, the child dimension, the parameter and the composition.
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
