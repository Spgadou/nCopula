addNode <- function(type,
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
    ## Class name
    name_short <- tolower(name_short)
    char1 <- strsplit(name_short, "")[[1]]
    char1[1] <- toupper(char1[1])
    char1 <- paste(char1, sep = "")
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
              contains = type)
    else
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
               contains = type)

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
        t@Laplace <- Laplace
        t@LaplaceInv <- LaplaceInv
        t@PGF <- PGF
        t@PGFInv <- PGFInv
      }
      else
      {
        t@Param <- "alpha"
        t@Laplace <- Laplace
        t@LaplaceInv <- LaplaceInv
      }
      t@simul <- simul
      t@theta <- vector("numeric")
      t@cop <- cop

      t
    })
    eval(parse(text = paste(".GlobalEnv$", toupper(name_short), " <- FF", sep = "")))
  }
}
