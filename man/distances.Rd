% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/CTfns.R, R/basicSCRfunctions.r, R/utilities.r
\name{distances}
\alias{distances}
\title{Calcluate distances.}
\usage{
distances(from, to)

distances(from, to)

distances(from, to)
}
\arguments{
\item{from}{A data frame or matrix with two columns, being the x- and y-coordinates of a set of points}

\item{to}{A data frame or matrix with two columns, being the x- and y-coordinates of a set of points}

\item{X, }{Y Matrices of coordinates.}

\item{X}{is an A by 2 matrix.}

\item{Y}{is a B by 2 matrix.}
}
\value{
Returns the matrix of distances.
}
\description{
Calculates euclidean distances between two sets of coordinates.

Returns all the distances between two sets of points, in a matix with as many rows as there are points 
 in \code{from} and as many columns as there are points in \code{t}..
}
\details{
Calculates distances between points in X and Y, whic are 
2-column matrices of coordinates (x- and y-coordinates in the two columns).
Returns a matrix with dimensions A=dim(X)[1] and B=dim(Y)[1], containing the distances.
}
\references{
This function written by Murray Efford
}
