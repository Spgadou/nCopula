% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/pCompCop.R
\name{pCompCop}
\alias{pCompCop}
\title{Distribution function of Mother class objects}
\usage{
pCompCop(structure, vector = FALSE, express = TRUE)
}
\arguments{
\item{structure}{object of class Mother.}

\item{vector}{logical. If false, returns a function or a character string with (u_1, u_2, ...) as arguments, else,
just (u).}

\item{express}{logical. If false, returns a function, else, a character string.}
}
\value{
The distribution function in the form of either a function or a character string.
}
\description{
Distribution function of a Mother class object.
}
\examples{
## Create the structure
structure <- LOG(0.5, NULL, list(GAMMA(1/30, c(5,6), NULL),
                              LOG(0.1, NULL, list(GAMMA(1/30, c(1,2), NULL),
                              GAMMA(1/30, c(3,4), NULL)))))

## Character string
pCompCop(structure, vector = TRUE, express = TRUE)
pCompCop(structure, vector = FALSE, express = TRUE)

## Function
pCompCop(structure, vector = TRUE, express = FALSE)
pCompCop(structure, vector = FALSE, express = FALSE)

}
