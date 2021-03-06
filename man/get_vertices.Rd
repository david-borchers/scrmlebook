% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/lcdplotting.R
\name{get_vertices}
\alias{get_vertices}
\title{Finds vertex index on graph made from raster}
\usage{
get_vertices(pts, rastermask)
}
\arguments{
\item{pts}{Matrix whose rows are (x,y) coordinates of points on raster}

\item{raster}{Raster object.}
}
\value{
Returns the row numbers of raster that correpond to pts. Note that pts must match exactly some 
coordinates of raster (use \code{closest_coords} to find closest coordinates if necessary).
}
\description{
Finds vertex index on graph made from raster, of points at coordinates pts. Vertex index is just the row of 
 the point in the raster object.
}
