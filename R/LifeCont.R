#' Contrat d'assurance vie discrète
#'
#' @description (...)
#' @param qx qx d'une table de mortalité
#' @export

assuranceVie <- function(qx, n = 1, x = 0, m = 0, delta = 0.04, g = 100000, kappa = 0.99){

  v <- exp(-delta)
  px <- 1 - qx

  FxBarre <- cumprod(px)
  FtBarre <- function(t)
  {
    if (x == 0)
      FxBarre[t+1]
    else
      FxBarre[t + x + 1] / FxBarre[x]
  }

  vk <- 1:n + m

  ## Valeurs possibles

  val.poss <- c(0, rev(g * v^vk)) * v^m

  dens.fun <- function(t) FtBarre(t) - FtBarre(t + 1)
  dens <- c(1 + FtBarre(m+n) - FtBarre(m), rev(dens.fun(vk - 1)))
  cdf <- cumsum(dens)

  res <- as.data.frame(cbind("z" = val.poss, "Densité" = dens, "CDF" = cdf))

  ## Mesures de risques

  VaRZ <- numeric(length(kappa))
  TVaRZ <- numeric(length(kappa))
  for (j in 1:length(kappa)){
  if (kappa[j] < FtBarre(n))
  {
    VaRZ[j] <- 0
    TVaRZ[j] <- (sum(val.poss[which(cdf > kappa[j])] * dens[which(cdf > kappa[j])])) / (1 - kappa[j])
  }

  else
  {
    VaRZ[j] <- val.poss[min(which(cdf > kappa[j]))]
    TVaRZ[j] <- (sum(val.poss[which(cdf > kappa[j])] * dens[which(cdf > kappa[j])]) + VaRZ[j] * (cdf[min(which(cdf > kappa[j])) - 1] - kappa[j])) / (1 - kappa[j])
  }}

  esp <- sum(val.poss * dens)
  var <- sum(val.poss^2 * dens) - esp^2
  sd <- sqrt(var)

  res.divers <- as.data.frame(cbind("Espérance" = esp, "SD" = sd, "Kappa" = kappa, "VaR" = VaRZ, "TVaR" = TVaRZ))


  list("Données" = res, "Mesures" = res.divers)
}


#' Contrat de rentes discrètes
#'
#' @description (...)
#' @param qx qx d'une table de mortalité
#' @export

rentes <- function(qx, n = 1, x = 0, m = 0,delta = 0.04, g = 10000, kappa = 0.99){

  v <- exp(-delta)
  px <- 1 - qx

  FxBarre <- cumprod(px)
  FtBarre <- function(t)
  {
    if (x == 0)
      FxBarre[t]
    else
      FxBarre[t + x] / FxBarre[x]
  }

  ## Valeurs possibles

  vk <- 0:(n-1) + m
  val.poss <- c(0, cumsum(g * v^vk))
  dens <- numeric(n)
  dens[n] <- FtBarre(n+m-1)
  for (i in 0:(n-2))
    dens[i+1] <- FtBarre(i+m) - FtBarre(i +m+ 1)

  dens <- c(1 - FtBarre(m), dens)

  cdf <- cumsum(dens)

  res <- as.data.frame(cbind("z" = val.poss, "Densité" = dens, "CDF" = cdf))

  ## Mesures

  esp <- sum(val.poss * dens)
  var <- sum(val.poss^2 * dens) - esp^2
  sd <- sqrt(var)

  VaRZ <- numeric(length(kappa))
  TVaRZ <- numeric(length(kappa))

  for (j in 1:length(kappa)){

  if (kappa[j] < 1 - FtBarre(m))
  {
    VaRZ[j] <- 0
    TVaRZ[j] <- sum(val.poss[which(cdf > kappa[j])] * dens[which(cdf > kappa[j])])
  }
  if (cdf[min(which(cdf > kappa[j]))] == 1)
  {
    VaRZ[j] <- val.poss[min(which(cdf > kappa[j]))]
    TVaRZ[j] <- VaRZ[j]
  }
  else
  {
    VaRZ[j] <- val.poss[min(which(cdf > kappa[j]))]
    TVaRZ[j] <- (sum(val.poss[which(cdf > kappa[j])] * dens[which(cdf > kappa[j])]) + VaRZ[j] * (cdf[min(which(cdf > kappa[j])) - 1] - kappa[j])) / (1 - kappa[j])
  }}

  res.divers <- as.data.frame(cbind("Espérance" = esp, "SD" = sd, "Kappa" = kappa, "VaR" = VaRZ, "TVaR" = TVaRZ))

  list("Données" = res, "Mesures" = res.divers)
}

#' SOA canadian mortality tables
#'
#' @description Table 2052.
#'
#' @usage data("t2052")
#'
#' @docType data
#'
#' @keywords datasets

"t2052"

#' SOA canadian mortality tables
#'
#' @description Table 2053.
#'
#' @usage data("t2053")
#'
#' @docType data
#'
#' @keywords datasets

"t2053"

#' SOA canadian mortality tables
#'
#' @description Table 2060.
#'
#' @usage data("t2060")
#'
#' @docType data
#'
#' @keywords datasets

"t2060"

#' SOA canadian mortality tables
#'
#' @description Table 2061.
#'
#' @usage data("t2061")
#'
#' @docType data
#'
#' @keywords datasets

"t2061"

#' SOA canadian mortality tables
#'
#' @description Table 2790.
#'
#' @usage data("t2790")
#'
#' @docType data
#'
#' @keywords datasets

"t2790"

#' SOA canadian mortality tables
#'
#' @description Table 2791.
#'
#' @usage data("t2791")
#'
#' @docType data
#'
#' @keywords datasets

"t2791"
