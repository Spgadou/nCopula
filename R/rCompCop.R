#' Density, Cdf, and Random Number Generator for Copulas Constructed Through Compounding
#'
#' @description The density, the CDF, and the random number generator for copulas construted
#' through compounding.
#' @param n Number of realisations
#' @param FUN Object of class Mother
#' @param level Number of imbrications
#'
#' @details rCompCop2 is more general (and easier to use) than rCompCop, but is slower.
#'
#' @author Simon-Pierre Gadoury
#'
#' @export

rCompCop <- compiler::cmpfun(function(n, FUN, level){

  if (FUN@type == "Child")
  {
    FUN@theta <- FUN@simul(n, FUN@parameter)
    laplace <- stringr::str_replace_all(FUN@Laplace, FUN@Param, FUN@parameter)
    th <- "-log(runif(FUN@dimension * n)) / FUN@theta"
    laplace <- parse(text = stringr::str_replace_all(laplace, "z", th))
    return(matrix(eval(laplace), ncol = FUN@dimension, nrow = n))
  }

  t <- FUN
  res <- list()
  FUN <- list()

  th.pos1 <- list()
  m.pos1 <- list()

  th.pos1[[1]] <- list()
  m.pos1[[1]] <- list()

  FUN[[1]] <- list()
  FUN[[1]][[1]] <- t


  for (i in 1:level)
  {
    FUN[[i + 1]] <- list()

    if (i == 1){

      FUN[[1]][[1]]@theta <- FUN[[1]][[1]]@simul(n, FUN[[1]][[1]]@parameter)

      FUN[[1]][[1]]@PGF <- stringr::str_replace_all(FUN[[1]][[1]]@PGF, FUN[[1]][[1]]@Param, FUN[[1]][[1]]@parameter)

      if (FUN[[1]][[1]]@type == "Mother"){
        typ <- numeric(FUN[[1]][[1]]@dimension)
        for (j in 1:length(FUN[[1]][[1]]@structure))
          typ[j] <- FUN[[1]][[1]]@structure[[j]]@type}

      th.pos1[[1]][[1]] <- which(typ == "Child")
      m.pos1[[1]][[1]] <- which(typ == "Mother")

      if (length(FUN[[1]][[1]]@arg) > 1 || {length(FUN[[1]][[1]]@arg) == 1 && FUN[[1]][[1]]@arg != 0})
      {
        for (j in 1:length(FUN[[1]][[1]]@arg))
        {
          laplace <- stringr::str_replace_all(FUN[[1]][[1]]@Laplace, FUN[[1]][[1]]@Param, FUN[[1]][[1]]@parameter)
          th <- "-log(runif(n)) / FUN[[1]][[1]]@theta"
          laplace <- parse(text = stringr::str_replace_all(laplace, "z", th))
          res[[FUN[[1]][[1]]@arg[j]]] <- eval(laplace)
        }
      }

      if (length(th.pos1[[1]][[1]]) != 0){
        for (j in 1:length(th.pos1[[1]][[1]]))
        {
          laplace <- stringr::str_replace_all(FUN[[1]][[1]]@structure[[th.pos1[[1]][[1]][j]]]@Laplace, FUN[[1]][[1]]@structure[[th.pos1[[1]][[1]][j]]]@Param, FUN[[1]][[1]]@structure[[th.pos1[[1]][[1]][j]]]@parameter)
          fbarre <- stringr::str_replace_all(FUN[[1]][[1]]@PGF, "z", laplace)

          argg <- FUN[[1]][[1]]@structure[[th.pos1[[1]][[1]][j]]]@arg

          th2 <- vapply(FUN[[1]][[1]]@theta, function(x) sum(FUN[[1]][[1]]@structure[[th.pos1[[1]][[1]][j]]]@simul(x, FUN[[1]][[1]]@structure[[th.pos1[[1]][[1]][j]]]@parameter)), rep(0, 1))

          res2 <- list()
          for (l in 1:length(argg))
          {
            th <- -log(runif(n)) / th2

            fbarre <- stringr::str_replace_all(fbarre, "z", "th")

            res[[argg[l]]] <- eval(parse(text = fbarre))
          }
        }}

      if (length(m.pos1[[1]][[1]]) != 0){
        for (j in 1:length(m.pos1[[1]][[1]]))
        {

          FUN[[1]][[1]]@structure[[m.pos1[[1]][[1]][j]]]@PGF <- stringr::str_replace_all(FUN[[1]][[1]]@structure[[m.pos1[[1]][[1]][j]]]@PGF, FUN[[1]][[1]]@structure[[m.pos1[[1]][[1]][j]]]@Param, FUN[[1]][[1]]@structure[[m.pos1[[1]][[1]][j]]]@parameter)

          FUN[[1]][[1]]@structure[[m.pos1[[1]][[1]][j]]]@PGF <- stringr::str_replace_all(FUN[[1]][[1]]@PGF, "z", FUN[[1]][[1]]@structure[[m.pos1[[1]][[1]][j]]]@PGF)

          FUN[[1]][[1]]@structure[[m.pos1[[1]][[1]][j]]]@theta <- vapply(FUN[[1]][[1]]@theta, function(x) sum(FUN[[1]][[1]]@structure[[m.pos1[[1]][[1]][j]]]@simul(x, FUN[[1]][[1]]@structure[[m.pos1[[1]][[1]][j]]]@parameter)), rep(0, 1))

          FUN[[i + 1]][[j]] <- FUN[[1]][[1]]@structure[[m.pos1[[1]][[1]][j]]]
        }}

    }

    else
    {
      th.pos1[[i]] <- list()
      m.pos1[[i]] <- list()

      for (k in 1:length(m.pos1[[i - 1]][[1]]))
      {

        typ <- numeric(FUN[[i]][[k]]@dimension)
        for (j in 1:length(FUN[[i]][[k]]@structure))
          typ[j] <- FUN[[i]][[k]]@structure[[j]]@type

        th.pos1[[i]][[k]] <- which(typ == "Child")
        m.pos1[[i]][[k]] <- which(typ == "Mother")

        if (FUN[[i]][[k]]@arg != 0)
        {
          for (j in 1:length(FUN[[i]][[k]]@arg))
          {
            laplace <- stringr::str_replace_all(FUN[[i]][[k]]@Laplace, FUN[[i]][[k]]@Param, FUN[[i]][[k]]@parameter)
            th <- "-log(runif(n)) / FUN[[i]][[k]]@theta"
            laplace <- parse(text = stringr::str_replace_all(laplace, "z", th))
            res[[FUN[[i]][[k]]@arg[j]]] <- eval(laplace)
          }
        }

        if (length(th.pos1[[i]][[k]]) != 0){
          for (j in 1:length(th.pos1[[i]][[k]]))
          {
            laplace <- stringr::str_replace_all(FUN[[i]][[k]]@structure[[th.pos1[[i]][[k]][j]]]@Laplace, FUN[[i]][[k]]@structure[[th.pos1[[i]][[k]][j]]]@Param, FUN[[i]][[k]]@structure[[th.pos1[[i]][[k]][j]]]@parameter)

            fbarre <- stringr::str_replace_all(FUN[[i]][[k]]@PGF, "z", laplace)

            argg <- FUN[[i]][[k]]@structure[[th.pos1[[i]][[1]][j]]]@arg

            th2 <- vapply(FUN[[i]][[k]]@theta, function(x) sum(FUN[[i]][[k]]@structure[[th.pos1[[i]][[k]][j]]]@simul(x, FUN[[i]][[k]]@structure[[th.pos1[[i]][[k]][j]]]@parameter)), rep(0, 1))

            res2 <- list()
            for (l in 1:length(argg))
            {
              th <- -log(runif(n)) / th2

              fbarre <- parse(text = stringr::str_replace_all(fbarre, "z", "th"))

              res[[argg[l]]] <- eval(fbarre)
            }

            #res[[argg[1]]] <- matrix(unlist(res2), ncol = FUN[[1]][[1]]@structure[[th.pos1[[1]][[1]][j]]]@dimension, nrow = n)
          }}

        if (length(m.pos1[[i]][[k]]) != 0){
          for (j in 1:length(m.pos1[[i]][[k]]))
          {

            FUN[[i]][[k]]@structure[[m.pos1[[i]][[k]][j]]]@PGF <- stringr::str_replace_all(FUN[[i]][[k]]@structure[[m.pos1[[i]][[k]][j]]]@PGF,
                                                                                           FUN[[i]][[k]]@structure[[m.pos1[[i]][[k]][j]]]@Param,
                                                                                           FUN[[i]][[k]]@structure[[m.pos1[[i]][[k]][j]]]@parameter)

            FUN[[i]][[k]]@structure[[m.pos1[[i]][[k]][j]]]@PGF <- stringr::str_replace_all(FUN[[i]][[k]]@PGF, "z", FUN[[i]][[k]]@structure[[m.pos1[[i]][[k]][j]]]@PGF)

            FUN[[i]][[k]]@structure[[m.pos1[[i]][[k]][j]]]@theta <- vapply(FUN[[i]][[k]]@theta, function(x) sum(FUN[[i]][[k]]@structure[[m.pos1[[i]][[k]][j]]]@simul(x, FUN[[i]][[k]]@structure[[m.pos1[[i]][[k]][j]]]@parameter)), rep(0, 1))

            FUN[[i + 1]][[j]] <- FUN[[i]][[k]]@structure[[m.pos1[[i]][[k]][j]]]

          }}
      }
    }

  }
  do.call(cbind, res)
})
