#' @slot Param The name of the parameter used
#' @slot parameter The value of the parameter
#' @slot dimension The dimension
#' @slot type The type of function (either child or mother)
#' @slot arguments The corresponding arguments (ex.: arguments 1 and 2 imply 'u1' and 'u2')
#' @slot structure The structure below the node of type 'Mother'
#' @slot Laplace Expression of the LST
#' @slot LaplaceInv Expression of the inverse LST
#' @slot PGF Expression of the pgf
#' @slot PGFInv Expression of the inverse pgf
#' @slot simul Fonction to sample from the distribution
#' @slot theta I don't know honestly
#' @slot cop Construct an Archimedean copula with this distribution
#' @slot Der Fonction to compute the expression of the 'k'th derivative of either the 'PGF', 'PGFInv', 'Laplace' or 'LaplaceInv'
#' @slot FUN Fonction to compute the function of the 'k'th derivative of either the 'PGF', 'PGFInv', 'Laplace' or 'LaplaceInv'
