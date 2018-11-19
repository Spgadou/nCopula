\name{AMH}
\alias{AMH}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{
%%  ~~function to do ... ~~
}
\description{
%%  ~~ A concise (1-5 lines) description of what the function does. ~~
}
\usage{
AMH(param, dim = 2L, density = FALSE)
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{param}{
%%     ~~Describe \code{param} here~~
}
  \item{dim}{
%%     ~~Describe \code{dim} here~~
}
  \item{density}{
%%     ~~Describe \code{density} here~~
}
}
\details{
%%  ~~ If necessary, more details than the description above ~~
}
\value{
%%  ~Describe the value returned
%%  If it is a LIST, use
%%  \item{comp1 }{Description of 'comp1'}
%%  \item{comp2 }{Description of 'comp2'}
%% ...
}
\references{
%% ~put references to the literature/web site here ~
}
\author{
%%  ~~who you are~~
}
\note{
%%  ~~further notes~~
}

%% ~Make other sections like Warning with \section{Warning }{....} ~

\seealso{
%% ~~objects to See Also as \code{\link{help}}, ~~~
}
\examples{
##---- Should be DIRECTLY executable !! ----
##-- ==>  Define data, use random,
##--	or do  help(data=index)  for the standard data sets.

## The function is currently defined as
function (param, dim = 2L, density = FALSE) 
{
    if (param < 0 || param >= 1) 
        stop("Wrong 'param' input")
    verif <- eval(parse(text = paste0(dim, "L")))
    if (!is.integer(verif) || dim <= 1) 
        stop("The dimension must be an integer greater than or equal to 2")
    phi <- "(alpha) / (exp(z) - (1 - alpha))"
    phi.inv <- "log((alpha + (z) * (1 - alpha)) / (z))"
    dep.param <- "alpha"
    param.th <- "alpha"
    th <- function(z, alpha) rgeom(z, alpha) + 1
    param <- as.character(param)
    if (density) {
        tt <- GEO(1/10, 1:dim, NULL)
        uu <- paste("u", 1:dim, sep = "")
        expr1 <- numeric(dim)
        for (i in 1:dim) expr1[i] <- stringr::str_replace_all(tt@Der("z", 
            1, "LaplaceInv"), "z", uu[i])
        expr1 <- paste("(", expr1, ")", sep = "", collapse = " * ")
        nu <- numeric(dim)
        for (i in 1:dim) nu[i] <- stringr::str_replace_all(tt@LaplaceInv, 
            "z", uu[i])
        nu <- paste("(", nu, ")", sep = "", collapse = " + ")
        expr2 <- stringr::str_replace_all(tt@Der("z", dim, "Laplace"), 
            "z", nu)
        densit <- paste("(", expr1, ") * (", expr2, ")", sep = "")
        densit <- stringr::str_replace_all(densit, "gamma", "alpha")
        new("amh", phi = phi, phi.inv = phi.inv, rBiv = rBiv, 
            theta = th, depend = dep.param, dimension = dim, 
            parameter = param, dens = densit, par.th = param.th, 
            name = "AMH copula")
    }
    else {
        new("amh", phi = phi, phi.inv = phi.inv, rBiv = rBiv, 
            theta = th, depend = dep.param, dimension = dim, 
            parameter = param, par.th = param.th, name = "AMH copula")
    }
  }
}
% Add one or more standard keywords, see file 'KEYWORDS' in the
% R documentation directory.
\keyword{ ~kwd1 }% use one of  RShowDoc("KEYWORDS")
\keyword{ ~kwd2 }% __ONLY ONE__ keyword per line
