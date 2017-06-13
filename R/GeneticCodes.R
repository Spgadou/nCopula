#' Obtain the genetic codes of a structure
#'
#' @param str The structure
#'
#' @return A list of of the structure's genetic codes.
#'
#' @author Simon-Pierre Gadoury
#'
#' @export

GeneticCodes <- function(str)
{
  e1 <- new.env(hash = TRUE, parent = parent.frame(), size = 10L)

  e1$ll <- list()
  e1$k <- 1
  e1$v <- list(c(0))

  FUN <- function(str, l = 1)
  {
    vk <- length(str@structure)
    type <- numeric(vk)
    for (i in 1:vk)
      type[i] <- str@structure[[i]]@type

    if (sum(str@arg) != 0)
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
          FUN(str@structure[[i]], l + 1)
        }
      }
    }
  }

  FUN(str)
  e1$ll
}
