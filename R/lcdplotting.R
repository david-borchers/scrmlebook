geommean = function(x) exp(mean(x)) # need to use link = list(noneuc="identity")

geomLCdist <- function (xy1, xy2, mask) {
  if (missing(xy1)) return('noneuc') # required by secr
  require(gdistance) # to load transition and geoCorrection functions
  if(is.element("noneuc",names(covariates(mask)))) {
    Sraster <- raster(mask, 'noneuc') # Make raster from mesh with covariate 'noneuc'
  } else if(is.element("noneuc.0",names(covariates(mask)))) {
    Sraster <- raster(mask, 'noneuc.0') # Make raster from mesh with covariate 'noneuc'
  } else stop("Got to have covariate named `noneuc` or `noneuc.0` on mask.")  
  # Calculate all the conductances, using mytransfun
  trans <- transition(Sraster, transitionFunction = geommean, directions = 16)
  # Adjust for distance difference between square and diagonal neighbours
  trans <- geoCorrection(trans)
  # calculate the least-cost distance (trans must be conductances)
  costDistance(trans, as.matrix(xy1), as.matrix(xy2))
}

arithmean = function(x) mean(exp(x)) # need to use (default) link = list(noneuc="log")

arithLCdist <- function (xy1, xy2, mask) {
  if (missing(xy1)) return('noneuc') # required by secr
  require(gdistance) # to load transition and geoCorrection functions
  if(is.element("noneuc",names(covariates(mask)))) {
    Sraster <- raster(mask, 'noneuc') # Make raster from mesh with covariate 'noneuc'
  } else if(is.element("noneuc.0",names(covariates(mask)))) {
    Sraster <- raster(mask, 'noneuc.0') # Make raster from mesh with covariate 'noneuc'
  } else stop("Got to have covariate named `noneuc` or `noneuc.0` on mask.")  
  # Calculate all the conductances, using mytransfun
  trans <- transition(Sraster, transitionFunction = arithmean, directions = 16)
  # Adjust for distance difference between square and diagonal neighbours
  trans <- geoCorrection(trans)
  # calculate the least-cost distance (trans must be conductances)
  costDistance(trans, as.matrix(xy1), as.matrix(xy2))
}

expmax = function(x) exp(max(x))

expmaxLCdist <- function (xy1, xy2, mask) {
  if (missing(xy1)) return('noneuc') # required by secr
  require(gdistance) # to load transition and geoCorrection functions
  if(is.element("noneuc",names(covariates(mask)))) {
    Sraster <- raster(mask, 'noneuc') # Make raster from mesh with covariate 'noneuc'
  } else if(is.element("noneuc.0",names(covariates(mask)))) {
    Sraster <- raster(mask, 'noneuc.0') # Make raster from mesh with covariate 'noneuc'
  } else stop("Got to have covariate named `noneuc` or `noneuc.0` on mask.")  
  # Calculate all the conductances, using mytransfun
  trans <- transition(Sraster, transitionFunction = expmax, directions = 16)
  # Adjust for distance difference between square and diagonal neighbours
  trans <- geoCorrection(trans)
  # calculate the least-cost distance (trans must be conductances)
  costDistance(trans, as.matrix(xy1), as.matrix(xy2))
}

expmaxcommdist <- function (xy1, xy2, mask) {
  if (missing(xy1)) return('noneuc') # required by secr
  require(gdistance) # to load transition and geoCorrection functions
  if(is.element("noneuc",names(covariates(mask)))) {
    Sraster <- raster(mask, 'noneuc') # Make raster from mesh with covariate 'noneuc'
  } else if(is.element("noneuc.0",names(covariates(mask)))) {
    Sraster <- raster(mask, 'noneuc.0') # Make raster from mesh with covariate 'noneuc'
  } else stop("Got to have covariate named `noneuc` or `noneuc.0` on mask.")  
  # Calculate all the conductances, using mytransfun
  trans <- transition(Sraster, transitionFunction = expmax, directions = 16)
  # Adjust for distance difference between square and diagonal neighbours
  trans <- geoCorrection(trans)
  # calculate the least-cost distance (trans must be conductances)
  commuteDistance(trans, as.matrix(xy1))
}

lcusageplot = function(fit,n=512,ac=NULL,mask=NULL,base="noneuc.0",lcdfun="geomLCdist",latent.class=1,maskcol=parula(40),
                       plotdets=FALSE,detpar=NULL,...) {
  require(pals)
  
  if(is.null(mask)) mask = fit$mask
  
  lcd.fun = match.fun(lcdfun)
  
  if(!is.element("noneuc.0",names(covariates(mask))))
    stop("Must have 'noneuc.0' as one of the mask covariates. It is not there.")
  if(!is.element(base,names(covariates(mask)))) {
    warning(paste("mask does not have a covariate called ",base,"; noneuc.0 being used instead."))
  }
  covariates(mask)$base = covariates(mask)$noneuc.0
  
  if(!is.null(ac)) {
    if(is.data.frame(ac)) ac = as.matrix(ac)
    if(!is.matrix(ac)) stop("ac must be a data frame or matrix")
    if(length(dim(ac))!=2) stop("ac must be a 2D data frame or matrix")
    n = dim(ac)[1]
  }
  
  for(i in 1:n){
    if(is.null(ac)) plotcovariate(mask,"base",col=maskcol,what="image")
    if(is.null(ac)) fromind = nearesttrap(unlist(locator(1)), mask)
    else fromind = nearesttrap(ac[i,], mask)
    from = c(x=mask$x[fromind],y=mask$y[fromind])
    from=matrix(from,ncol=2)
    dists = lcd.fun(from,mask,mask)
    dfn <- secr:::getdfn(fit$detectfn)
    detpars = detectpar(fit,byclass=TRUE)
    if(is.list(detpars)) { # this means there is more than one latent class
      nclass = length(detpars)
      if(nclass<latent.class) {
        warning(paste("latent.class argument > number of classes, latent.class=",nclass," used",sep=""))
        latent.class = nclass
      }
      detpars = unlist(detpars[[latent.class]])
    } else {
      detpars = unlist(detpars)
    }
    p = dfn(dists[1,],detpars)
    covariates(mask)$p = p/sum(p)
    plotcovariate(mask,"p",what="image",...)
    points(from,col="white",pch=19)
    points(from,cex=0.8,pch=19)
    if(is.null(ac)) waitclick = unlist(locator(1))
    
    if(plotdets) {
      cams = traps(fit$capthist)
      plot(cams,add=TRUE,detpar=detpar)
    }
  }
  
  invisible(mask)
}


#' @title Least-cost path plot
#'
#' @description
#'  Plots 
#'  Requires packages secr, raster, gdistance, fields, igraph
#'  
#' @param mask A mask with covariate `noneuc` (obtained using \code{predictDsurface})
#' @param transitionFunction A user-defined function for calculating conductance between points 
#' (see function \code{transition} of the \code{gdistance} package).
#' 'userdfn1' in the \code{secr} vignette 'secr-noneuclidean.pdf', for example.)
#' @param covariate The name (a character) of a covariate to plot as an image plot
#' (typically because 'noneuc' depends on it).
#' @param n number of least-cost paths to plot
#' @param linecol color of least cost path line
#' @param directions parameter of function \code{transition}
#' @param symm parameter of function \code{transition}
#' @param directed parameter of function \code{transition}
#' @param ... arguments to function \code{plotcovariate}, which is used to plot the image.
#' 
#' @examples 
#' lcpathplot(Dhat.Nemegt,"stdBC",userdfun=function(x) 1/x,n=5,contour=FALSE,col=c("orange","brown"),key=FALSE, linecol="white",xaxt="n",yaxt="n",xlab="",ylab="",bty="n",asp=1)
#' @export 
lcpathplot = function(mask,transitionFunction,type="noneuc",n=512,froms=NULL,tos=NULL,linecol="white", 
                      directions=16, symm=TRUE,directed=FALSE,lwd=1,...) {
  
  if(!is.element("noneuc.0",names(covariates(mask))))
    stop("Must have 'noneuc.0' as one of the mask covariates. It is not there.")
  if(!is.element(type,c("noneuc","sigma")))
    stop(paste("Invalid type: '",type,"'"," It must be `noneuc` or `sigma`.",sep=""))
  rastermask = raster(mask,"noneuc.0") # make raster with covariates(mask)$noneuc.0 as values of pixels
  
  transfun=match.fun(transitionFunction)
  
  coords = coordinates(rastermask) # lookup table for vertex coordinates
  # secr models conductance (sigma) with a log link. Line below assumes that $noneuc.0 is on linear predictor scale 
  # tr1<-transition(rastermask,transitionFunction=function(x) exp(lp(x)),directions=directions,symm=symm)
  tr1<-transition(rastermask,transitionFunction=transfun,directions=directions,symm=symm)
  tr1CorrC<-geoCorrection(tr1,type="c",multpl=FALSE,scl=FALSE)
  
  if(type=="noneuc") plotcovariate(mask,"noneuc.0",...)
  else {
    covariates(mask)$sigma = exp(covariates(mask)$noneuc.0)
    plotcovariate(mask,"sigma",...)
  }
  
  if(!is.null(froms)) {
    if(is.data.frame(froms)) froms = as.matrix(froms)
    if(!is.matrix(froms)) stop("froms must be a data frame or matrix")
    if(length(dim(froms))!=2) stop("froms must be a 2D data frame or matrix")
    n = dim(froms)[1]
    if(is.null(tos)) stop("If you specify froms, you must also specify tos")
    if(is.data.frame(tos)) tos = as.matrix(tos)
    if(!is.matrix(tos)) stop("tos must be a data frame or matrix")
    if(length(dim(tos))!=2) stop("tos must be a 2D data frame or matrix")
    if(dim(froms)[1]!=n) stop("dim(tos) must be same as dim(froms)" )
  }
  
  dists = rep(NA,n) # to keep least-cost path distances in
  if(is.null(froms)) froms = tos = data.frame(x=rep(NA,n), y=rep(NA,n)) # to keep coords in
  for(i in 1:n) {
    if(is.null(froms)) fromind = nearesttrap(unlist(locator(1)), mask)
    else fromind = nearesttrap(froms[i,], mask)
    if(is.null(tos)) toind = nearesttrap(unlist(locator(1)), mask)
    else toind = nearesttrap(tos[i,], mask)
    from = c(x=mask$x[fromind],y=mask$y[fromind])
    to = c(x=mask$x[toind],y=mask$y[toind])
    froms[i,] = from
    tos[i,] = to
    from=matrix(from,ncol=2)
    to=matrix(to,ncol=2)
    #    npts=dim(from)[1]
    #    nptsto=dim(to)[1]
    #    if(nptsto != npts) stop("Must have same number of points in 'from' as in 'to'.")
    #    if(npts>1) pts = closest_coords(from[i,],to[i,],rastermask)
    #    else pts = closest_coords(from,to,rastermask)
    pts = closest_coords(from,to,rastermask)
    vpts = get_vertices(pts,rastermask)
    
    trmat=summary(transitionMatrix(tr1CorrC))
    #cbind(trmat,1/trmat$x)
    #    rel=data.frame(from=trmat$i,to=trmat$j,weight=1/trmat$x)
    rel=data.frame(from=trmat$i,to=trmat$j,weight=trmat$x)
    #rel
    g = graph_from_data_frame(rel,directed=directed,vertices=NULL)
    #    attributes(g)$noneuc.0=1/trmat$x
    #    E(g)$weight=1/trmat$x
    attributes(g)$noneuc.0=trmat$x
    E(g)$weight=trmat$x
    #vertices = as_data_frame(g, what="vertices")
    #edges = as_data_frame(g, what="edges")
    svert=which(names(V(g))==vpts[1])
    evert=which(names(V(g))==vpts[2])
    # NB: Need to invert E(g) so that higher values lead to shorter distances:
    spath=as.numeric(names(shortest_paths(g,from=svert,to=evert,weights=1/E(g)$weight)$vpath[[1]]))
    dists[i]=igraph:::distances(g,v=svert,to=evert,weights=attributes(g)$noneuc.0)
    
    nppts=length(spath)
    segments(coords[spath[-nppts],1],coords[spath[-nppts],2],coords[spath[-1],1],coords[spath[-1],2],col=linecol,lwd=lwd)
    points(coords[spath[c(1,nppts)],],pch=19,col="white",cex=1.5)
    points(coords[spath[c(1,nppts)],],pch=19,col=c("green","red"),cex=0.75)
  }
  
  invisible(list(from=froms,to=tos,dists=dists))
}


#' @title Finds closest coordinates on raster to two points
#'
#' @description
#'  Uses function over() from package sp to overlay points on raster and return closest raster coordinates
#'  
#' @param from pair of coordinates (x,y) from which to start
#' @param to pair of coordinates (x,y) to get to
#' @param rastermask Raster object (typically created from mask by something like 
#' rastermask = raster(mask,"noneuc"))
#' 
#' @return Returns the coordinates of the closest point on the raster, as a matrix with two columns (x,y), 
#' named s1 and s2, with first row corresponding to 'from' coordinates, and second row corresponding to 'to' 
#' coordinates.
#' @export closest_coords
#' 
closest_coords=function(from,to,rastermask){
  ends=SpatialPoints(rbind(from,to))
  grid=as(rastermask, 'SpatialGrid') 
  xy=over(ends,grid)
  return(coordinates(grid)[xy,])
}


#' @title Finds vertex index on graph made from raster
#'
#' @description
#'  Finds vertex index on graph made from raster, of points at coordinates pts. Vertex index is just the row of 
#'  the point in the raster object.
#'  
#' @param pts Matrix whose rows are (x,y) coordinates of points on raster
#' @param raster Raster object.
#' 
#' @return Returns the row numbers of raster that correpond to pts. Note that pts must match exactly some 
#' coordinates of raster (use \code{closest_coords} to find closest coordinates if necessary).
#' 
#' @export get_vertices
#' 
get_vertices = function(pts,rastermask){
#  target = nearest.raster.point(pts[,1],pts[,2],as.im(rastermask),indices=FALSE)
  target = nearest.raster.point(pts[,1],pts[,2],as.im(as.matrix(rastermask)),indices=FALSE)
  coords = coordinates(rastermask) # lookup table from index produced by transition() to coordinates
  npts = dim(pts)[1]
  vert = rep(NA,npts)
  for(i in 1:npts){
    #    vert[i] = which(coords[,1]==target$x[i] & coords[,2]==target$y[i])
    dst = sqrt((coords[,1]-target$x[i])^2 + (coords[,2]-target$y[i])^2)
    vert[i] = which(dst == min(dst))[1]
  }
  return(vert)
}

#' @title Creates the igraph of a mask object
#'
#' @description
#'  Creates an igraph object with a vertex for each mask point and edges to neighbours, weighted according 
#'  to the function \code{cfun}, using the mask covariate \code{cname}.
#'  
#'  Requires packages raster, gdistance, igraph
#'  
#' @param mask \code{secr} mask object. Must have covariate called 'noneuc' containing cost
#' @param cname Name of variable to use in cost calculation
#' @param cfun Cost function name. Defauls to "mean"
#' @param directed If TRUE, use directed graph for transition between cells, else undirected 
#' @param symm If TRUE, cost is same in both directions 
#' @param conductance If TRUE, cfun must calculate conductance (the inverse of cost), not cost.
#' 
#' @examples 
#' ny=4; nx=4 # dimensions of mask
#' # set costs (NA is "hole" - nothing there & can't go there):
#' costs=c(100,100,100,100,NA,100,100,NA,1,NA,100,1,1,1,1,1) 
#' rmesh=data.frame(x=rep(1:nx,ny),y=rep(1:ny,rep(nx,ny)),noneuc=costs) # make data frame with coords and costs
#' 
#' rmask=read.mask(data=rmesh,columns="noneuc")  # make mask with covariate 'noneuc' containing cost
#' ig=make_igraph(rmask,"noneuc")
#' plot(ig, edge.label=round(E(g)$weight, 3))
#' 
#' cfun=function(x) exp(diff(x)) # asymmetric cost function
#' 
#' ig=make_igraph(rmask,"noneuc",cfun=cfun)
#' plot(ig, edge.label=round(E(g)$weight, 3))
#' 
#' ig=make_igraph(rmask,"noneuc",cfun=cfun,directed=TRUE)
#' plot(ig, edge.label=round(E(g)$weight, 3))
#' 
#' ig=make_igraph(rmask,"noneuc",cfun=cfun,directed=TRUE)
#' plot(ig, edge.label=round(E(g)$weight, 3))
#' 
#' ig=make_igraph(rmask,"noneuc",cfun=cfun,directed=TRUE,symm=FALSE)
#' plot(ig, edge.label=round(E(g)$weight, 3))
#' 
#' @export make_igraph
#' 
make_igraph = function(mask,cname,cfun="mean",directions=8,directed=FALSE,symm=TRUE,conductance=FALSE) {
  require(raster)
  require(spatstat)
  
  
  if(!is.element(cname,names(covariates(mask))))
    stop(paste("'",cname,"'"," is not the name of one of the mask covariates.",sep=""))
  if(mask$y[1]<mask$y[dim(mask)[1]]) mask$y = max(mask$y)-mask$y+1 # reorder so count down from top (y-axis increasing downwards)
  rastermask = raster(mask,cname) # make raster with covariates(mask)$cname as values of pixels
  
  f=match.fun(cfun)
  #tr1<-transition(rastermask,transitionFunction=function(x) 1/mean(x),directions=8)
  #tr1<-transition(rastermask,transitionFunction=function(x) 1/exp(diff(x)),directions=8,symm=FALSE)
  if(conductance) tr1<-transition(rastermask,transitionFunction=function(x) f(x),directions=directions,symm=symm)
  else tr1<-transition(rastermask,transitionFunction=function(x) 1/f(x),directions=directions,symm=symm)
  tr1CorrC<-geoCorrection(tr1,type="c",multpl=FALSE,scl=FALSE)
  #costs1<-costDistance(tr1CorrC,pts)
  
  from = as.data.frame(mask)
  to = as.data.frame(mask)
  pts = closest_coords(from,to,rastermask)
  vpts = get_vertices(pts,rastermask)
  
  trmat=summary(tr1CorrC)
  if(conductance) rel=data.frame(from=trmat$i,to=trmat$j,weight=trmat$x)
  else rel=data.frame(from=trmat$i,to=trmat$j,weight=1/trmat$x)
  if(directed) g = graph_from_data_frame(rel,directed=TRUE,vertices=NULL)
  else g = graph_from_data_frame(rel,directed=FALSE,vertices=NULL)
  attributes(g)$noneuc=1/trmat$x
  E(g)$weight=1/trmat$x
  
  return(g)  
}



#' @title Least-cost distance summary
#'
#' @description
#'  Tabulates least-cost distance between smallest and largest `noneuc` values on mask
#'  Requires packages secr, raster, gdistance, fields, igraph
#'  
#' @param masknu A mask with covariate `noneuc` (obtained using \code{predictDsurface})
#' @param userfn A user-defined function for calculating least-cost distances (see function 
#' 'userdfn1' in the \code{secr} vignette 'secr-noneuclidean.pdf', for example.)
#' @param covariate The name (a character) of a covariate to summarise with the distances 
#' (typically because 'noneuc' depends on it).
#' 
noneucdists = function(masknu,userfn,covariate=NULL) {
  if(!is.element("noneuc",names(covariates(masknu)))) stop("Need a mask covariate called `noneuc'.")
  if(!is.character(covariate)) stop("Argument `covariate' must be a character.")
  if(!is.null(covariate) & !is.element(covariate,names(covariates(masknu)))) 
    stop(paste("No covariate called",covariate,"in the mask"))
  f = match.fun(userfn)
  uNU = range(covariates(masknu)$noneuc)
  loNUind = which(covariates(masknu)$noneuc==uNU[1])[1]
  hiNUind = which(covariates(masknu)$noneuc==uNU[2])[1]
  lonu = covariates(masknu)$noneuc[loNUind]
  hinu = covariates(masknu)$noneuc[hiNUind]
  if(!is.null(covariate)) {
    covcol = which(names(covariates(masknu))==covariate)
    locov = covariates(masknu)[,covcol][loNUind]
    hicov = covariates(masknu)[,covcol][hiNUind]
  }
  # make mask with two points, distance 2 apart
  dets = trap.builder(frame=data.frame(x=c(-2,2),y=c(0,0)),method="all")
  msk = make.mask(dets,buffer=2, spacing=1)
  # calculate dists for lo-lo, lo-hi, hi-lo, hi-hi
  covariates(msk) = data.frame(noneuc=rep(c(lonu,hinu,hinu,lonu,lonu,rep(0,dim(msk)[1]-5))))
  alldists = f(as.matrix(msk),as.matrix(msk),msk)
  dists = c(alldists[1,2],alldists[2,3],alldists[3,4],alldists[4,5])
  distnames = c("Lo.noneuc-Hi.noneuc","Hi.noneuc-Hi.noneuc","Hi.noneuc-Lo.noneuc","Lo.noneuc-Lo.noneuc")
  if(!is.null(covariate)) {
    dists=matrix(rep(dists,5),nrow=5)
    dists[2,] = c(lonu,hinu,hinu,lonu)
    dists[3,] = c(hinu,hinu,lonu,lonu)
    dists[4,] = c(locov,hicov,hicov,locov)
    dists[5,] = c(hicov,hicov,locov,locov)
    row.names(dists) = c("Distance","from.noneuc","to.noneuc",paste("from.",covariate,sep=""),paste("to.",covariate,sep=""))
    colnames(dists) = distnames
  } else {
    dists=matrix(rep(dists,3),nrow=3)
    dists[2,] = c(lonu,hinu,hinu,lonu)
    dists[3,] = c(hinu,hinu,lonu,lonu)
    row.names(dists) = c("Distance","from.noneuc","to.noneuc")
    colnames(dists) = distnames
  }
  dists
}


#' @title Add leas cost path to plot
#'
#' @description
#'  Calculates and plots the least-cost path between two points according 
#'  to the function \code{cfun}, using the mask covariate \code{cname}.
#'  
#'  Requires packages raster, gdistance, igraph
#'  
#' @param mask \code{secr} mask object. Must have covariate called 'noneuc' containing cost
#' @param cname Name of variable to use in cost calculation
#' @param cfun Cost function name. Defauls to "mean"
#' @param directed If TRUE, use directed graph for transition between cells, else undirected 
#' @param symm If TRUE, cost is same in both directions 
#' @param conductance If TRUE, cfun must calculate conductance (the inverse of cost), not cost.
#' 
#' @details Invisibly returns an object of class \code{SpatialLines} defining the least cost path.
#' 
#' @export add_lcpath
#' 
add_lcpath = function(mask,cname,cfun="mean",directed=FALSE,symm=TRUE,conductance=FALSE,...) {
  require(raster)
  require(spatstat)
  require(gdistance)
  
  
  if(!is.element(cname,names(covariates(mask))))
    stop(paste("'",cname,"'"," is not the name of one of the mask covariates.",sep=""))
  rastermask = raster(mask,cname) # make raster with covariates(mask)$cname as values of pixels
  
  f=match.fun(cfun)
  #tr1<-transition(rastermask,transitionFunction=function(x) 1/mean(x),directions=8)
  #tr1<-transition(rastermask,transitionFunction=function(x) 1/exp(diff(x)),directions=8,symm=FALSE)
  if(conductance) tr1<-transition(rastermask,transitionFunction=function(x) f(x),directions=8,symm=symm)
  else tr1<-transition(rastermask,transitionFunction=function(x) 1/f(x),directions=8,symm=symm)
  lcpath = shortestPath(tr1,from,to,output="SpatialLines")
  lines(lcpath,...)

invisible(lcpath)

}

