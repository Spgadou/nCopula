#' Obtain the Genetic Codes of a Structure
#'
#' @description Function to obtain the list of all genetic codes of a structure.
#'
#' @param structure an object of class Mother (the structure)
#'
#' @return A list of the structure's genetic codes.
#'
#' @examples
#' ## Create the structure
#' structure <- GEO(0.5, NULL, list(GAMMA(1/30, c(5,6), NULL),
#'                            GEO(0.1, NULL, list(GAMMA(1/30, c(1,2), NULL),
#'                                                GAMMA(1/30, c(3,4), NULL)))))
#' ## Get the genetic codes
#' GeneticCodes(structure)
#'
#' @author Simon-Pierre Gadoury
#'
#' @export

GeneticCodes <- function(structure)
{
  e1 <- new.env(hash = TRUE, parent = parent.frame(), size = 10L)

  e1$ll <- list()
  e1$k <- 1
  e1$v <- list(c(0))

  FUN <- function(structure, l = 1)
  {
    vk <- length(structure@structure)
    type <- numeric(vk)
    for (i in 1:vk)
      type[i] <- structure@structure[[i]]@type

    if (sum(structure@arg) != 0)
    {
      e1$ll[[e1$k]] <- c(e1$v[[l]], 0)
      e1$k <- e1$k + 1
    }

    if (sum(type == "Mother") == 0)
    {
      for (i in 1:vk)
      {
        e1$ll[[e1$k]] <- c(e1$v[[l]], i)
        e1$k <- e1$k + 1
      }
      e1$ll
    }
    else
    {
      for (i in 1:vk)
      {
        if (type[i] == "Child")
        {
          e1$ll[[e1$k]] <- c(e1$v[[l]], i)
          e1$k <- e1$k + 1
        }
        else
        {
          e1$v[[l + 1]] <- c(e1$v[[l]], i)
          FUN(structure@structure[[i]], l + 1)
        }
      }
    }
  }

  FUN(structure)
  e1$ll
}
