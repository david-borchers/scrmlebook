% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/lcdplotting.R
\name{make_igraph}
\alias{make_igraph}
\title{Creates the igraph of a mask object}
\usage{
make_igraph(mask, cname, cfun = "mean", directions = 8,
  directed = FALSE, symm = TRUE, conductance = FALSE)
}
\arguments{
\item{mask}{\code{secr} mask object. Must have covariate called 'noneuc' containing cost}

\item{cname}{Name of variable to use in cost calculation}

\item{cfun}{Cost function name. Defauls to "mean"}

\item{directed}{If TRUE, use directed graph for transition between cells, else undirected}

\item{symm}{If TRUE, cost is same in both directions}

\item{conductance}{If TRUE, cfun must calculate conductance (the inverse of cost), not cost.}
}
\description{
Creates an igraph object with a vertex for each mask point and edges to neighbours, weighted according 
 to the function \code{cfun}, using the mask covariate \code{cname}.
 
 Requires packages raster, gdistance, igraph
}
\examples{
ny=4; nx=4 # dimensions of mask
# set costs (NA is "hole" - nothing there & can't go there):
costs=c(100,100,100,100,NA,100,100,NA,1,NA,100,1,1,1,1,1) 
rmesh=data.frame(x=rep(1:nx,ny),y=rep(1:ny,rep(nx,ny)),noneuc=costs) # make data frame with coords and costs

rmask=read.mask(data=rmesh,columns="noneuc")  # make mask with covariate 'noneuc' containing cost
ig=make_igraph(rmask,"noneuc")
plot(ig, edge.label=round(E(g)$weight, 3))

cfun=function(x) exp(diff(x)) # asymmetric cost function

ig=make_igraph(rmask,"noneuc",cfun=cfun)
plot(ig, edge.label=round(E(g)$weight, 3))

ig=make_igraph(rmask,"noneuc",cfun=cfun,directed=TRUE)
plot(ig, edge.label=round(E(g)$weight, 3))

ig=make_igraph(rmask,"noneuc",cfun=cfun,directed=TRUE)
plot(ig, edge.label=round(E(g)$weight, 3))

ig=make_igraph(rmask,"noneuc",cfun=cfun,directed=TRUE,symm=FALSE)
plot(ig, edge.label=round(E(g)$weight, 3))

}
