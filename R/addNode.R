#' Add a distribution
#'
#' @param type the type of the node (either 'Mother', 'Child' or 'Both')
#' @param pp the parameter of the distribution used in the character strings
#' @param name_short the short name of the distribution (ex.: 'log' for the logarithmic distribution)
#' @param name_long the long name of the distribution (ex.: 'logarithmic' for the logarithmic distribution)
#' @param Laplace the LST of the distribution (character), where
#' @param LaplaceInv the inverse LST of the distribution (character)
#' @param PGF the pgf of the distribution (character), if type is 'Mother' or 'Both'
#' @param PGFInv the inverse pgf of the distribution (character), if type is 'Mother' or 'Both'
#' @param simul function to sample from the distribution
#' @param cop function to create a corresponding Archimedean copula class (ex.: for a GEO, it is an AMH copula), can be NULL
#'
#' @author Simon-Pierre Gadoury
#'
#' @examples
#' addNode(type = "Child",
#' pp = "gamma",
#' name_short = "pois",
#' name_long = "poisson",
#' Laplace = "exp(gamma * (exp(-(z)) - 1))",
#' LaplaceInv = "-log(log(z) / gamma + 1)",
#' NULL,
#' NULL,
#' simul = function(n, gamma) rpois(n, gamma),
#' NULL)
#'
#' POIS(5, 1:2, NULL)
#'
#' @export

addNode <- function(type,
                    pp,
                    name_short,
                    name_long,
                    Laplace,
                    LaplaceInv,
                    PGF = NULL,
                    PGFInv = NULL,
                    simul,
                    cop = NULL)
{
  if (type != "Mother" && type != "Child")
    stop("The type should be either 'Child' or 'Mother'")
  else
  {

    name_long <- tolower(name_long)
    char1 <- strsplit(name_long, "")[[1]]
    char1[1] <- toupper(char1[1])
    name_long <- paste(char1, collapse = "")

    ## Class name
    name_short <- tolower(name_short)
    char1 <- strsplit(name_short, "")[[1]]
    char1[1] <- toupper(char1[1])
    char1 <- paste(char1, collapse = "")
    name_class <- paste(char1, "_", type, sep = "")

    ## Create the class
    if (type == "Child")
      setClass(name_class,
             list(name = "character",
                  type = "character",
                  dimension = "numeric",
                  parameter = "numeric",
                  arg = "numeric",
                  obj = "character",
                  Param = "character",
                  Laplace = "character",
                  LaplaceInv = "character",
                  PGF = "character",
                  PGFInv = "character",
                  simul = "function",
                  theta = "numeric",
                  LTheta = "character",
                  cop = "function"),
              contains = type, where = .GlobalEnv)
    else if (type == "Mother")
      setClass(name_class,
               list(name = "character",
                    type = "character",
                    dimension = "numeric",
                    parameter = "numeric",
                    arg = "numeric",
                    structure = "list",
                    obj = "character",
                    Param = "character",
                    Laplace = "character",
                    LaplaceInv = "character",
                    PGF = "character",
                    PGFInv = "character",
                    simul = "function",
                    theta = "numeric",
                    PM = "character",
                    cop = "function"),
               contains = type, where = .GlobalEnv)
    else if (type == "Both")
    {
      type <- "Mother"
      name_short <- tolower(name_short)
      char1 <- strsplit(name_short, "")[[1]]
      char1[1] <- toupper(char1[1])
      char1 <- paste(char1, collapse = "")
      name_class <- paste(char1, "_", type, sep = "")

      setClass(name_class,
               list(name = "character",
                    type = "character",
                    dimension = "numeric",
                    parameter = "numeric",
                    arg = "numeric",
                    structure = "list",
                    obj = "character",
                    Param = "character",
                    Laplace = "character",
                    LaplaceInv = "character",
                    PGF = "character",
                    PGFInv = "character",
                    simul = "function",
                    theta = "numeric",
                    PM = "character",
                    cop = "function"),
               contains = type, where = .GlobalEnv)

      type <- "Child"
      name_short <- tolower(name_short)
      char1 <- strsplit(name_short, "")[[1]]
      char1[1] <- toupper(char1[1])
      char1 <- paste(char1, collapse = "")
      name_class <- paste(char1, "_", type, sep = "")
      setClass(name_class,
               list(name = "character",
                    type = "character",
                    dimension = "numeric",
                    parameter = "numeric",
                    arg = "numeric",
                    obj = "character",
                    Param = "character",
                    Laplace = "character",
                    LaplaceInv = "character",
                    PGF = "character",
                    PGFInv = "character",
                    simul = "function",
                    theta = "numeric",
                    LTheta = "character",
                    cop = "function"),
               contains = type, where = .GlobalEnv)
    }

    FF <- compiler::cmpfun(function(par, unif, struc)
    {
      if (length(unique(unif)) != length(unif))
        stop("The 'unif' argument must be composed of different values")

      if (is.null(struc))
      {
        t <- new(name_class, parameter = par, arg = unif, dimension = length(unif), name = paste(name_long, " distribution", sep = ""), type = type, obj = name_short)
      }

      else
      {
        if (class(struc) != "list")
          stop("The argument 'struc' must be a list")

        if (is.null(unif))
          t <- new(name_class, parameter = par, dimension = length(struc), structure = struc, arg = 0, paste(name_long, " distribution", sep = ""), type = type, obj = name_short)
        else
          t <- new(name_class, parameter = par, dimension = length(struc) + length(unif), structure = struc, arg = unif, paste(name_long, " distribution", sep = ""), type = type, obj = name_short)
      }

      if (t@type == "Mother")
      {
        t@Param <- "gamma"
        t@Laplace <- stringr::str_replace_all(Laplace, pp, "gamma")
        t@LaplaceInv <- stringr::str_replace_all(LaplaceInv, pp, "gamma")
        t@PGF <- stringr::str_replace_all(PGF, pp, "gamma")
        t@PGFInv <- stringr::str_replace_all(PGFInv, pp, "gamma")
      }
      else
      {
        t@Param <- "alpha"
        t@Laplace <- stringr::str_replace_all(Laplace, pp, "alpha")
        t@LaplaceInv <- stringr::str_replace_all(LaplaceInv, pp, "alpha")
      }
      t@simul <- simul
      t@theta <- vector("numeric")

      if (is.null(cop))
        t@cop <- function() "No known Archimedean copula related to this distribution"
      else
        t@cop <- cop

      t
    })
    eval(parse(text = paste(".GlobalEnv$", toupper(name_short), " <- FF", sep = "")))
  }
}
